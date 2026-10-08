"""Temporal-holdout validation of the deployment hotspots (Test 8).

    run()                               -> tables, results/TEMPORAL_HOLDOUT_RESULTS.md, draft figure, processed table
    make_view(T, ...)                   -> a dated input view (every input dated after T removed)
    score_view(view, condition, measurement, level) -> the unchanged engine's frame on that view

Pre-specified in docs/ANALYSIS_PLAN_TEMPORAL_HOLDOUT.md (written before any holdout outcome was computed).

Question: do regions the engine would have flagged with data available at a past cutoff T turn out to be where Long
COVID / ME/CFS research activity (trials, new trial sites, NIH awards) appeared after T? That is a PROXY outcome:
later research activity supports the engine's research-readiness side; it is not evidence of clinical value,
diagnostic yield or patient benefit. The desert side is built to point at places with little activity, so for it the
checkable statement is gap persistence (flagged deserts stayed under-served), reported separately.

The scoring engine is not forked. `opportunity.make_combo` / `combo_frame` and `metric_link.reach_for` run unchanged on
a *view* in which `opportunity.spatial()` (facility / provider base) and `opportunity._features()` (burden and desert
inputs) are replaced by date-filtered copies (context manager `engine_view`). `make_view(None)` reproduces the stored
`deployment_opportunities` composites (tested). Leave-one-state-out normalisation swaps only the percentile-rank
functions (`reference_normalisation`), so the same formula runs with a different normalisation reference.

Everything is ecological: places, registries and facilities, never people.
"""
from __future__ import annotations

import json
import math
import sys
import time
import warnings
from contextlib import contextmanager
from dataclasses import dataclass, field, replace
from functools import lru_cache

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import rankdata

from .. import config
from ..config import FIGURES, TABLES, UNKNOWN, utc_now_iso
from ..facilities.registry import SEP
from ..geography import features as F
from ..provenance import add_provenance
from ..store import read_table, write_table
from . import metric_link as ML
from . import opportunity as O

PRODUCER = "measure_it.scoring.temporal_holdout"
TABLE = "temporal_holdout_rankings"
PLAN = "docs/ANALYSIS_PLAN_TEMPORAL_HOLDOUT.md"
RESULTS_MD = config.RESULTS / "TEMPORAL_HOLDOUT_RESULTS.md"
PRIMARY_T = 2021
CUTOFFS = (2021, 2020, 2019)
QUERY = O.PRIMARY["condition"]                 # long_covid_or_me_cfs
MEASUREMENT = O.PRIMARY["measurement"]         # wearable_autonomic_activity_monitoring
QUERY_CONDS = ("long_covid", "me_cfs")
CMS_CODE = "51"
CMS_LATEST_FINAL = 2022
N_BOOT = 1000
N_PERM = 1000
KS = (10, 25, 100)
RADIUS_KM = 50
# variants (plan section 3): A = clean (Long COVID burden undefined before T); B = Long COVID burden held at its
# post-period value (NOT a clean holdout); static = A with providers / facilities NOT date-filtered (leakage bound)
VARIANTS = {
    "A_clean": {"lc_burden": "undefined", "static_context": False, "clean": True,
                "label": "clean: Long COVID burden undefined before T; ME/CFS CMS proxy of year T"},
    "B_lc_post_burden": {"lc_burden": "post", "static_context": False, "clean": False,
                         "label": "NOT a clean holdout: Long COVID burden held at its current (post-period) HPS value"},
    "A_static_context": {"lc_burden": "undefined", "static_context": True, "clean": False,
                         "label": "leakage bound: variant A with providers and facilities NOT date-filtered (today's "
                                  "NPPES / facility snapshot)"},
}
COMPONENT_SCORES = ["burden_pct", "vulnerability_pct", "diagnostic_desert_pct", "clinic_capacity_pct",
                    "research_readiness_pct"]
# pre-specified expected direction of each score for the research-activity outcomes (plan section 5.3)
EXPECTED = {
    "research_readiness_pct": "positive", "clinic_capacity_pct": "positive", "burden_pct": "none",
    "vulnerability_pct": "negative", "diagnostic_desert_pct": "negative",
    "composite_equal": "weakly positive", "composite_evidence_weighted": "weakly positive",
    "composite_study_partner": "positive", "composite_capacity_led": "positive",
    "composite_access_gap": "around or below 0.5 (desert ranking)", "composite_burden_led": "none",
    "composite_equity_led": "around or below 0.5 (desert ranking)", "composite_burden_only": "none",
    "composite_saturation_adjusted": "weakly positive",
    "baseline_population_within_50km": "positive", "baseline_population_total": "positive",
    "baseline_prior_trials_50km": "positive", "baseline_random": "0.5",
}
SCORE_KIND = {**{c: "component" for c in COMPONENT_SCORES}}
GUARDRAIL = ("Proxy outcome: research activity that appeared after the cutoff is a research-readiness signal, not "
             "evidence of clinical value, diagnostic yield or patient benefit. Ecological (places, registries, "
             "facilities); nothing describes a person. A ranked region is a candidate deployment opportunity.")


def log(msg: str) -> None:
    print(f"[temporal_holdout] {msg}", file=sys.stderr, flush=True)


def cutoff_ts(T: int) -> pd.Timestamp:
    return pd.Timestamp(f"{int(T)}-12-31")


# ------------------------------------------------------------------------------------------------------------------
# dated inputs
# ------------------------------------------------------------------------------------------------------------------

@lru_cache(maxsize=1)
def trial_dates() -> pd.DataFrame:
    """nct_id, first_post (study_first_post_date), start (start_date; month precision -> first of month)."""
    ct = read_table("clinical_trials", columns=["nct_id", "study_first_post_date", "start_date"])
    return pd.DataFrame({"nct_id": ct["nct_id"].astype(str),
                         "first_post": pd.to_datetime(ct["study_first_post_date"], errors="coerce", format="mixed"),
                         "start": pd.to_datetime(ct["start_date"], errors="coerce", format="mixed")}).drop_duplicates("nct_id")


def trials_before(T: int | None, rule: str = "first_post") -> set:
    """NCT ids whose `rule` date is on or before T-12-31 (all trials when T is None). A trial without the date is
    never pre-period (it could not be dated)."""
    d = trial_dates()
    if T is None:
        return set(d["nct_id"])
    return set(d.loc[d[rule].notna() & (d[rule] <= cutoff_ts(T)), "nct_id"])


def trials_after(T: int, rule: str = "first_post") -> set:
    d = trial_dates()
    return set(d.loc[d[rule].notna() & (d[rule] > cutoff_ts(T)), "nct_id"])


@lru_cache(maxsize=1)
def relevant_trials() -> pd.DataFrame:
    """Literal-view trial x condition links (trial_conditions.condition_literal_match), deduplicated."""
    tc = read_table("trial_conditions", columns=["nct_id", "condition_id", "condition_literal_match"])
    tc = tc[tc["condition_literal_match"].astype(bool)].drop_duplicates(["nct_id", "condition_id"])
    return tc[["nct_id", "condition_id"]].assign(nct_id=lambda d: d["nct_id"].astype(str)).reset_index(drop=True)


@lru_cache(maxsize=1)
def npi_enumeration() -> pd.Series:
    p = read_table("providers", columns=["npi", "enumeration_date"])
    return pd.Series(pd.to_datetime(p["enumeration_date"]).to_numpy(), index=p["npi"].astype(str))


@lru_cache(maxsize=1)
def facility_dates() -> pd.DataFrame:
    """Earliest dated source per facility: NPI enumeration, HRSA site added to scope, first-posted date of any trial at
    the facility, first NIH fiscal year at the facility (plan section 2)."""
    f = read_table("facilities", columns=["facility_id", "npi_object_ids", "hrsa_site_object_ids"])
    out = pd.DataFrame({"facility_id": f["facility_id"].astype(str)}).set_index("facility_id")
    npis = f[["facility_id", "npi_object_ids"]].assign(
        npi=f["npi_object_ids"].fillna("").astype(str).str.split(SEP)).explode("npi")
    npis = npis[npis["npi"].fillna("").ne("")]
    npis["npi"] = npis["npi"].str.replace("npi:", "", regex=False)
    npis["d"] = npis["npi"].map(npi_enumeration())
    out["npi_first"] = npis.groupby("facility_id")["d"].min()
    h = read_table("facilities__hrsa", columns=["site_bphc_number", "site_added_to_scope_date"])
    hd = pd.Series(pd.to_datetime(h["site_added_to_scope_date"]).to_numpy(), index="hrsa_site:" + h["site_bphc_number"].astype(str))
    hs = f[["facility_id", "hrsa_site_object_ids"]].assign(
        s=f["hrsa_site_object_ids"].fillna("").astype(str).str.split(SEP)).explode("s")
    hs = hs[hs["s"].fillna("").ne("")]
    hs["d"] = hs["s"].map(hd)
    out["hrsa_first"] = hs.groupby("facility_id")["d"].min()
    ft = read_table("facility_trials", columns=["facility_id", "nct_id"])
    ft["d"] = ft["nct_id"].astype(str).map(trial_dates().set_index("nct_id")["first_post"])
    out["trial_first"] = ft.groupby("facility_id")["d"].min()
    fn = read_table("facility_nih_projects", columns=["facility_id", "fiscal_year"])
    out["nih_first_fy"] = fn.groupby("facility_id")["fiscal_year"].min().astype(float)
    return out


def facility_exists(T: int, facility_ids) -> np.ndarray:
    fd = facility_dates().reindex(pd.Index(facility_ids).astype(str))
    c = cutoff_ts(T)
    ok = (fd["npi_first"] <= c) | (fd["hrsa_first"] <= c) | (fd["trial_first"] <= c) | (fd["nih_first_fy"] <= T)
    return ok.fillna(False).to_numpy(bool)


@lru_cache(maxsize=1)
def _spatial_npi_dates() -> np.ndarray:
    """Enumeration date per individual NPI in the order of opportunity.spatial() (same read and filter)."""
    p = read_table("providers", columns=["npi", "entity_type", "state", "county_fips", "lat", "lon",
                                         "specialty_groups", "enumeration_date"])
    p = p[(p["entity_type"].astype(str) == "1") & p["lat"].notna() & p["specialty_groups"].fillna("").ne("")]
    return pd.to_datetime(p["enumeration_date"]).to_numpy()


@lru_cache(maxsize=1)
def _npi_state_groups_dated() -> pd.DataFrame:
    """opportunity._npi_state_groups() with each NPI's enumeration date (same filter)."""
    p = read_table("providers", columns=["npi", "entity_type", "state_fips", "specialty_groups", "enumeration_date"])
    p = p[(p["entity_type"].astype(str) == "1") & p["state_fips"].notna() & p["specialty_groups"].fillna("").ne("")]
    x = p[["npi", "state_fips", "enumeration_date"]].assign(
        group=p["specialty_groups"].astype(str).str.split(SEP)).explode("group")
    return x[x["group"].fillna("").ne("")].reset_index(drop=True)


@lru_cache(maxsize=1)
def cms_burden() -> pd.DataFrame:
    """CMS MMD condition 51 (ME/CFS level-C proxy), unsmoothed actual, all ages, every claims year, 2024 geographies."""
    b = read_table("geo_condition_burden", columns=["burden_row_id", "geo_id", "geo_level", "condition_id", "measure_id",
                                                    "measure_label", "value", "value_unit", "year", "period",
                                                    "source_name", "geo_match_status"])
    b = b[(b["condition_id"] == "me_cfs") & (b["measure_id"] == f"cms_mmd_{CMS_CODE}_prev_actual_all")
          & (b["geo_match_status"] == "matched_2024")]
    return b.drop_duplicates(["geo_level", "geo_id", "year"]).reset_index(drop=True)


# ------------------------------------------------------------------------------------------------------------------
# the dated view and the engine swap
# ------------------------------------------------------------------------------------------------------------------

@dataclass
class View:
    T: int | None
    variant: str
    S: O.Spatial
    npi_groups: pd.DataFrame
    raw: dict = field(default_factory=dict)       # (cond, level) -> frame of raw desert inputs aligned to universe
    features: pd.DataFrame | None = None
    stats: dict = field(default_factory=dict)


@contextmanager
def engine_view(view: View):
    """Run the unchanged engine on `view`: opportunity.spatial / _features / _npi_state_groups swapped, input caches
    cleared on entry and exit."""
    saved = (O.spatial, O._features, O._npi_state_groups, O.capacity_basis_for)
    O._item_incidence.cache_clear()
    ML._zcta_covered.cache_clear()
    O.spatial = lambda: view.S
    if view.features is not None:
        O._features = lambda: view.features
    O._npi_state_groups = lambda: view.npi_groups
    if view.T is not None:
        # activity capacity (Medicare 2024 billing, clinic_capacity_basis since 2026-10-07) post-dates every cutoff:
        # a dated view keeps the implementer-density capacity built from its own dated facility base
        O.capacity_basis_for = lambda measurement_id: "implementer_density"
    try:
        yield view
    finally:
        O.spatial, O._features, O._npi_state_groups, O.capacity_basis_for = saved
        O._item_incidence.cache_clear()
        ML._zcta_covered.cache_clear()


def ref_percentile(x, ref: np.ndarray) -> np.ndarray:
    """Percentile rank against a reference subset (plan section 7): reference rows keep the engine definition
    (average rank among the non-null reference / its size); a non-reference row gets the percentile it would have if
    added to the reference (average rank among reference + itself / (size + 1)). NaN stays NaN."""
    x = np.asarray(x, dtype=float)
    ok = ~np.isnan(x)
    r = np.sort(x[ok & ref])
    n = len(r)
    out = np.full(x.shape, np.nan)
    if n == 0:
        return out
    lt = np.searchsorted(r, x, side="left").astype(float)
    le = np.searchsorted(r, x, side="right").astype(float)
    eq = le - lt
    inref = ref & ok
    out[inref] = (lt[inref] + (eq[inref] + 1) / 2) / n
    nr = ~ref & ok
    out[nr] = (lt[nr] + (eq[nr] + 2) / 2) / (n + 1)
    return out


@contextmanager
def reference_normalisation(ref: np.ndarray):
    """Swap the engine's (and the desert's) percentile-rank functions for the reference version on county-length
    arrays (len(ref)); every other array is ranked by the original functions."""
    n = len(ref)
    orig_pr, orig_rows, orig_f = O.percentile_rank, O.pct_rank_rows, F.percentile_rank

    def pr(x):
        a = np.asarray(x, dtype=float)
        return ref_percentile(a, ref) if a.ndim == 1 and len(a) == n else orig_pr(x)

    def rows(X):
        X = np.atleast_2d(np.asarray(X, dtype=float))
        if X.shape[1] != n:
            return orig_rows(X)
        return np.vstack([ref_percentile(X[i], ref) for i in range(X.shape[0])])

    def fpr(s):
        if len(s) == n:
            return pd.Series(ref_percentile(s.to_numpy(float), ref), index=s.index)
        return orig_f(s)

    O.percentile_rank, O.pct_rank_rows, F.percentile_rank = pr, rows, fpr
    try:
        yield
    finally:
        O.percentile_rank, O.pct_rank_rows, F.percentile_rank = orig_pr, orig_rows, orig_f


@contextmanager
def _filtered_reads(filters: dict):
    """geography.features counters read their inputs through `read_table`; filter named tables for one call."""
    orig = F.read_table

    def rt(name, *a, **k):
        d = orig(name, *a, **k)
        return filters[name](d) if name in filters else d

    F.read_table = rt
    try:
        yield
    finally:
        F.read_table = orig


def county_points() -> pd.DataFrame:
    return O.spatial().levels["county"].geo[["geo_id", "lat", "lon"]].reset_index(drop=True)


def trial_pool_counts(ncts: set, conds=QUERY_CONDS, label: str = "union") -> tuple[pd.DataFrame, pd.DataFrame]:
    """Distinct relevant trials (literal view, conditions `conds`, NCT ids in `ncts`) with a U.S. site in the county or
    within 50 / 100 km of its internal point (county) and in the state (state): geography.features.trial_counts
    unchanged, on a filtered trial_conditions (conditions relabelled `label` so the union counts once)."""
    keep = set(map(str, ncts))

    def flt(d):
        d = d[d["condition_id"].isin(list(conds)) & d["nct_id"].astype(str).isin(keep)].copy()
        d["condition_id"] = label
        return d.drop_duplicates(["nct_id", "condition_id"]) if len(d) else d

    counties = county_points()
    with _filtered_reads({"trial_conditions": flt}):
        try:
            c = F.trial_counts(counties)
        except ValueError:   # no trial at all (pd.concat of nothing)
            c = F.Counts(pd.DataFrame(), pd.DataFrame())
    out = []
    for lvl in (c.county, c.state):
        if len(lvl) and label in lvl.index.get_level_values(0):
            out.append(lvl.xs(label, level=0))
        else:
            out.append(pd.DataFrame())
    return out[0], out[1]


@lru_cache(maxsize=64)
def pool_counts(T: int | None, rule: str, conds: tuple, label: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Cached trial_pool_counts: rule 'pre' = first posted <= T; 'first_post' / 'start' = after T by that date."""
    ncts = trials_before(T) if rule == "pre" else trials_after(T, rule)
    return trial_pool_counts(ncts, conds=conds, label=label)


def build_raw(view: View) -> None:
    """Raw desert inputs per (condition, level) for the view: burden value, relevant providers per 100k, relevant
    trials in pool. T=None: the stored geo_condition_features inputs (and burden)."""
    f0 = O._features()
    for level in O.LEVELS:
        geo = view.S.levels[level].geo
        for cid in QUERY_CONDS:
            base = f0[(f0["condition_id"] == cid) & (f0["geo_level"] == level)].set_index("geo_id").reindex(geo["geo_id"])
            raw = pd.DataFrame({"geo_id": geo["geo_id"].to_numpy()})
            if view.T is None:
                raw["burden_value"] = base["burden_value"].to_numpy(float)
                raw["relevant_providers_per_100k"] = base["relevant_providers_per_100k"].to_numpy(float)
                tcol = "trials_in_geo_or_within_50km_n" if level == "county" else "trials_in_geo_n"
                raw["trials_pool_n"] = base[tcol].to_numpy(float)
            else:
                # trials of THIS condition first posted <= T (the desert's trial term is per condition)
                ctc, cts = pool_counts(view.T, "pre", (cid,), cid)
                src = ctc if level == "county" else cts
                tcol = "trials_in_geo_or_within_50km_n" if level == "county" else "trials_in_geo_n"
                raw["trials_pool_n"] = (src[tcol].reindex(geo["geo_id"]).fillna(0).to_numpy(float)
                                        if len(src) and tcol in src else np.zeros(len(geo)))
                with engine_view(view):
                    _, rate = O.relevant_provider_counts(cid, level)
                raw["relevant_providers_per_100k"] = rate
                if cid == "me_cfs":
                    cb = cms_burden()
                    cb = cb[(cb["geo_level"] == level) & (cb["year"] == float(view.T))].set_index("geo_id")
                    raw["burden_value"] = cb["value"].reindex(geo["geo_id"]).to_numpy(float)
                    raw["burden_row_id"] = cb["burden_row_id"].reindex(geo["geo_id"]).to_numpy()
                    raw["burden_measure_label"] = cb["measure_label"].reindex(geo["geo_id"]).to_numpy()
                    raw["burden_period"] = cb["period"].reindex(geo["geo_id"]).to_numpy()
                elif VARIANTS[view.variant]["lc_burden"] == "post":
                    raw["burden_value"] = base["burden_value"].to_numpy(float)
                else:
                    raw["burden_value"] = np.nan
            view.raw[(cid, level)] = raw


def build_features(view: View) -> pd.DataFrame:
    """geo_condition_features rows for the query's members, with the view's burden and desert (desert percentiles by
    geography.features._desert_components, the geography module's own function)."""
    f0 = O._features()
    frames = []
    for level in O.LEVELS:
        geo = view.S.levels[level].geo
        for cid in QUERY_CONDS:
            f = f0[(f0["condition_id"] == cid) & (f0["geo_level"] == level)].set_index("geo_id").reindex(geo["geo_id"])
            f = f.reset_index()
            raw = view.raw[(cid, level)]
            if view.T is not None:
                lc_undefined = cid == "long_covid" and VARIANTS[view.variant]["lc_burden"] == "undefined"
                f["burden_value"] = raw["burden_value"].to_numpy(float)
                has = f["burden_value"].notna()
                if cid == "me_cfs":
                    f["burden_row_id"] = raw["burden_row_id"].to_numpy()
                    f["burden_measure_label"] = raw["burden_measure_label"].to_numpy()
                    f["burden_period"] = raw["burden_period"].to_numpy()
                    f["burden_year"] = float(view.T)
                    f["burden_measure_id"] = f"cms_mmd_{CMS_CODE}_prev_actual_all"
                    f["burden_evidence_level"] = np.where(has, "C", "D")
                    f["burden_defined_level"] = "C"
                    f["burden_unknown_reason"] = np.where(has, None, f"no CMS MMD condition {CMS_CODE} value for this "
                                                                     f"geography in claims year {view.T}")
                elif lc_undefined:
                    for c in ("burden_row_id", "burden_measure_label", "burden_period", "burden_measure_id",
                              "burden_source_resolution"):
                        f[c] = None
                    f["burden_year"] = np.nan
                    f["burden_ci_low"] = np.nan
                    f["burden_ci_high"] = np.nan
                    f["burden_evidence_level"] = "D"
                    f["burden_defined_level"] = "D"
                    f["burden_inherited"] = False
                    f["burden_unknown_reason"] = (f"no Long COVID burden measure existed on or before {view.T} "
                                                  "(Household Pulse Survey estimates start in 2022)")
                f["relevant_providers_per_100k"] = raw["relevant_providers_per_100k"].to_numpy(float)
                tcol = "trials_in_geo_or_within_50km_n" if level == "county" else "trials_in_geo_n"
                f[tcol] = raw["trials_pool_n"].to_numpy(float)
            g = pd.DataFrame({"burden_value": raw["burden_value"].to_numpy(float),
                              "svi_overall": f["svi_overall"].to_numpy(float),
                              "relevant_providers_per_100k": raw["relevant_providers_per_100k"].to_numpy(float),
                              "trials": raw["trials_pool_n"].to_numpy(float), "in_analysis_universe": True})
            comp = F._desert_components(g, "burden_value", "relevant_providers_per_100k", "trials")
            f["desert_pct_burden"] = comp["b"].to_numpy()
            f["desert_pct_vulnerability"] = comp["v"].to_numpy()
            f["desert_pct_low_providers"] = comp["p"].to_numpy()
            f["desert_pct_low_trials"] = comp["t"].to_numpy()
            f["diagnostic_desert"] = comp[["b", "v", "p", "t"]].mean(axis=1, skipna=False).to_numpy()
            f["diagnostic_desert_access_only"] = comp[["v", "p", "t"]].mean(axis=1, skipna=False).to_numpy()
            if view.T is not None:
                f["desert_burden_resolution"] = np.where(f["burden_value"].notna(), "county" if level == "county"
                                                         else "state", None)
                if cid == "long_covid" and VARIANTS[view.variant]["lc_burden"] == "post" and level == "county":
                    f["desert_burden_resolution"] = "state (inherited)"
            f["in_analysis_universe"] = True
            frames.append(f)
    return pd.concat(frames, ignore_index=True)


def make_view(T: int | None, variant: str = "A_clean") -> View:
    """Dated input view (plan section 2). T=None: the current inputs (desert rebuilt from the stored raw inputs)."""
    S0 = O.spatial()
    if T is None:
        v = View(T=None, variant="current", S=S0, npi_groups=O._npi_state_groups())
    else:
        static = VARIANTS[variant]["static_context"]
        b0 = S0.base
        pre = trials_before(T)
        ft = b0.ft[b0.ft["nct_id"].astype(str).isin(pre)].copy()
        fn = b0.fn[b0.fn["fiscal_year"].astype(float) <= T].copy()
        G = b0.G.copy()
        exists = facility_exists(T, b0.fac["facility_id"].astype(str))
        if not static:
            G[~exists] = False
        base = replace(b0, G=G, ft=ft, fn=fn)
        dates = _spatial_npi_dates()
        assert len(dates) == len(S0.npi_bits), "provider order differs from opportunity.spatial()"
        npi_ok = ~(dates > np.datetime64(cutoff_ts(T)))
        bits = S0.npi_bits if static else np.where(npi_ok, S0.npi_bits, 0)
        ng = _npi_state_groups_dated()
        if not static:
            ng = ng[ng["enumeration_date"] <= cutoff_ts(T)]
        S = replace(S0, base=base, npi_bits=bits)
        v = View(T=T, variant=variant, S=S, npi_groups=ng[["npi", "state_fips", "group"]].reset_index(drop=True))
        v.stats.update({
            "facilities_clinic_candidates": int(len(b0.fac)), "facilities_existing_at_T": int(exists.sum()),
            "facilities_without_pre_T_date": int((~exists).sum()),
            "facility_trial_rows": int(len(b0.ft)), "facility_trial_rows_pre_T": int(len(ft)),
            "facility_nih_rows": int(len(b0.fn)), "facility_nih_rows_pre_T": int(len(fn)),
            "individual_npis": int(len(dates)), "individual_npis_enumerated_after_T": int((~npi_ok).sum()),
            "providers_date_filtered": not static, "facilities_date_filtered": not static,
        })
    build_raw(v)
    v.features = build_features(v)
    return v


def score_view(view: View, condition: str = QUERY, measurement: str = MEASUREMENT, level: str = "county",
               ref: np.ndarray | None = None) -> pd.DataFrame:
    """The unchanged engine (make_combo -> combo_frame, no Monte Carlo) on the view; `ref` = LOSO reference mask."""
    with engine_view(view):
        if ref is None:
            cb = O.make_combo(condition, measurement, level)
            return O.combo_frame(cb, mc=False)
        with reference_normalisation(ref):
            view.features = build_features(view)
            cb = O.make_combo(condition, measurement, level)
            ev = O.evaluate(cb)
            comps = O.named_composites(cb, ev)
        view.features = build_features(view)   # restore the ordinary normalisation for later calls
        out = pd.DataFrame({"geo_id": cb.geo["geo_id"].to_numpy(), "state_fips": cb.geo["state_fips"].to_numpy(),
                            "rank_eligible": cb.eligible})
        for k, c in comps.items():
            out[f"composite_{k}"] = c
        out["rank_equal"] = O.rank_rows(comps["equal"], cb.eligible)[0]
        return out


# ------------------------------------------------------------------------------------------------------------------
# post-period outcomes
# ------------------------------------------------------------------------------------------------------------------

@lru_cache(maxsize=1)
def nih_first_fy() -> pd.Series:
    pj = read_table("nih_projects", columns=["core_project_num", "fiscal_year"])
    return pj.groupby("core_project_num")["fiscal_year"].min().astype(float)


def nih_new_counts(T: int, conds=QUERY_CONDS) -> tuple[pd.Series, pd.Series]:
    """Relevant (precision view) NIH core projects first funded after T, per county / state of the awardee
    (geography.features.nih_counts, unchanged, on filtered links)."""
    pj = read_table("nih_projects", columns=["appl_id", "core_project_num"])
    first = pj["core_project_num"].map(nih_first_fy())
    new_appl = set(pj.loc[first > T, "appl_id"])

    def flt(d):
        d = d[d["condition_id"].isin(list(conds)) & d["appl_id"].isin(new_appl)].copy()
        d["condition_id"] = "union"
        return d

    with _filtered_reads({"nih_project_conditions": flt}):
        c = F.nih_counts()
    get = lambda x: x.xs("union", level=0)["nih_core_projects_n"] if len(x) and "union" in x.index.get_level_values(0) \
        else pd.Series(dtype=float)  # noqa: E731
    return get(c.county), get(c.state)


@lru_cache(maxsize=1)
def _facility_points() -> pd.DataFrame:
    f = read_table("facilities", columns=["facility_id", "lat", "lon", "county_fips", "state_fips"])
    return f.assign(facility_id=f["facility_id"].astype(str)).set_index("facility_id")


def new_facility_counts(T: int, conds=QUERY_CONDS) -> tuple[pd.Series, pd.Series, dict]:
    """Facilities with >= 1 relevant (literal view) trial first posted after T and none on or before T, counted in
    the county or within 50 km of its internal point (county) and in the state (state)."""
    ft = read_table("facility_trials", columns=["facility_id", "nct_id", "condition_ids_literal"])
    ft = ft[O.M._has(ft["condition_ids_literal"], list(conds))]
    ft["facility_id"] = ft["facility_id"].astype(str)
    pre = trials_before(T)
    post = trials_after(T)
    had_pre = set(ft.loc[ft["nct_id"].astype(str).isin(pre), "facility_id"])
    new = sorted(set(ft.loc[ft["nct_id"].astype(str).isin(post), "facility_id"]) - had_pre)
    pts = _facility_points().reindex(new)
    geo = O.spatial().levels["county"].geo
    ok = pts["lat"].notna().to_numpy()
    a = F.radius_matrix(geo["lat"].to_numpy(), geo["lon"].to_numpy(), pts["lat"].to_numpy()[ok],
                        pts["lon"].to_numpy()[ok], RADIUS_KM)
    gi = pd.Series(np.arange(len(geo)), index=geo["geo_id"])
    mem_r = pts["county_fips"].map(gi).to_numpy()[ok]
    from scipy import sparse
    has = ~pd.isna(mem_r)
    member = sparse.csr_matrix((np.ones(int(has.sum())), (mem_r[has].astype(int), np.flatnonzero(has))),
                               shape=a.shape)
    cnt = np.asarray(((a + member) > 0).sum(axis=1)).ravel().astype(float)
    county = pd.Series(cnt, index=geo["geo_id"].to_numpy())
    state = pts["state_fips"].value_counts()
    return county, state, {"new_facilities": len(new), "new_facilities_geocoded": int(ok.sum()),
                           "facilities_with_pre_T_relevant_trial": len(had_pre)}


def outcomes(T: int) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """County and state post-period outcomes (plan section 4)."""
    geo = O.spatial().levels["county"].geo
    sgeo = O.spatial().levels["state"].geo
    cid = geo["geo_id"]
    oc = pd.DataFrame({"geo_id": cid.to_numpy()})
    os_ = pd.DataFrame({"geo_id": sgeo["geo_id"].to_numpy(), "state_fips": sgeo["state_fips"].to_numpy()})
    post = trials_after(T)
    info = {"post_trials_any_condition": len(post)}
    rel = relevant_trials()
    for key, conds, rule in (("o1", QUERY_CONDS, "first_post"), ("o1start", QUERY_CONDS, "start"),
                             ("o2", ("me_cfs",), "first_post"), ("prior", QUERY_CONDS, "pre")):
        ncts = trials_before(T) if rule == "pre" else trials_after(T, rule)
        c, s = pool_counts(T, rule, tuple(conds), "union")
        oc[f"{key}_trials_pool_n"] = c["trials_in_geo_or_within_50km_n"].reindex(cid).fillna(0).to_numpy() \
            if len(c) else 0.0
        oc[f"{key}_trials_in_county_n"] = c["trials_in_geo_n"].reindex(cid).fillna(0).to_numpy() \
            if len(c) and "trials_in_geo_n" in c else 0.0
        os_[f"{key}_trials_in_state_n"] = s["trials_in_geo_n"].reindex(os_["state_fips"]).fillna(0).to_numpy() \
            if len(s) else 0.0
        info[f"{key}_relevant_trials"] = int(rel[rel["condition_id"].isin(list(conds))
                                                 & rel["nct_id"].isin(ncts)]["nct_id"].nunique())
    nc, ns = nih_new_counts(T)
    oc["o3_nih_new_core_in_county_n"] = nc.reindex(cid).fillna(0).to_numpy()
    os_["o3_nih_new_core_in_state_n"] = ns.reindex(os_["state_fips"]).fillna(0).to_numpy()
    fc, fs, finfo = new_facility_counts(T)
    oc["o4_new_facilities_pool_n"] = fc.reindex(cid).fillna(0).to_numpy()
    os_["o4_new_facilities_in_state_n"] = fs.reindex(os_["state_fips"]).fillna(0).to_numpy()
    info.update(finfo)
    info["o3_nih_new_core_projects"] = int(nc.sum())
    # O5 burden input persistence
    cb = cms_burden()
    cc = cb[cb["geo_level"] == "county"].pivot_table(index="geo_id", columns="year", values="value")
    oc["o5_cms51_T"] = cc.get(float(T), pd.Series(dtype=float)).reindex(cid).to_numpy()
    oc["o5_cms51_latest_final"] = cc.get(float(CMS_LATEST_FINAL), pd.Series(dtype=float)).reindex(cid).to_numpy()
    oc["population_total"] = geo["population_total"].to_numpy(float)
    oc["population_within_50km"] = geo["population_within_50km"].to_numpy(float)
    oc["state_fips"] = geo["state_fips"].to_numpy()
    os_["population_total"] = sgeo["population_total"].to_numpy(float)
    return oc, os_, info


OUTCOMES = {  # outcome key -> (binary column rule, count column, description)
    "O1": ("o1_trials_pool_n", "any relevant trial (Long COVID or ME/CFS) first posted after T with a U.S. site in the "
                               "county or within 50 km (PRIMARY)"),
    "O1_in_county": ("o1_trials_in_county_n", "as O1, site inside the county"),
    "O1_start": ("o1start_trials_pool_n", "as O1, trials partitioned by start date"),
    "O2": ("o2_trials_pool_n", "any ME/CFS trial first posted after T in the county or within 50 km"),
    "O3": ("o3_nih_new_core_in_county_n", "any relevant NIH core project first funded after T at an awardee in the "
                                          "county"),
    "O4": ("o4_new_facilities_pool_n", "any facility with a relevant post-T trial and no relevant pre-T trial in the "
                                       "county or within 50 km"),
}


# ------------------------------------------------------------------------------------------------------------------
# metrics (vectorised over bootstrap resamples / permutations)
# ------------------------------------------------------------------------------------------------------------------

def auroc_from_ranks(r: np.ndarray, y: np.ndarray) -> np.ndarray:
    """AUROC per row: r (B x n) average ranks of the score, y (B x n) binary."""
    y = y.astype(float)
    n1 = y.sum(axis=1)
    n0 = y.shape[1] - n1
    with np.errstate(invalid="ignore", divide="ignore"):
        a = ((r * y).sum(axis=1) - n1 * (n1 + 1) / 2) / (n1 * n0)
    a[(n1 == 0) | (n0 == 0)] = np.nan
    return a


def auroc(score, y) -> float:
    s = np.asarray(score, float)
    return float(auroc_from_ranks(rankdata(s)[None, :], np.asarray(y)[None, :])[0])


def _order_positions(score: np.ndarray, geo_ids: np.ndarray) -> np.ndarray:
    """Position of each row in the ranking (score descending, ties by geo_id), 0 = top."""
    order = np.lexsort((geo_ids, -score))
    pos = np.empty(len(score), dtype=np.int64)
    pos[order] = np.arange(len(score))
    return pos


def _topk_rate(pos_s: np.ndarray, y_s: np.ndarray, k: int) -> np.ndarray:
    """Share positive among the k best-positioned rows of each resample (rows of pos_s / y_s)."""
    k = min(k, pos_s.shape[1])
    idx = np.argpartition(pos_s, k - 1, axis=1)[:, :k]
    return np.take_along_axis(y_s, idx, axis=1).mean(axis=1)


def _spearman_rows(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    ra, rb = rankdata(a, axis=1), rankdata(b, axis=1)
    ra = ra - ra.mean(axis=1, keepdims=True)
    rb = rb - rb.mean(axis=1, keepdims=True)
    with np.errstate(invalid="ignore", divide="ignore"):
        return (ra * rb).sum(axis=1) / np.sqrt((ra ** 2).sum(axis=1) * (rb ** 2).sum(axis=1))


def _ci(x: np.ndarray) -> tuple[float, float]:
    x = np.asarray(x, float)
    x = x[~np.isnan(x)]
    if len(x) < 10:
        return np.nan, np.nan
    return float(np.percentile(x, 2.5)), float(np.percentile(x, 97.5))


def point_metrics(score: np.ndarray, y: np.ndarray, count: np.ndarray, geo_ids: np.ndarray) -> dict:
    pos = _order_positions(score, geo_ids)
    n = len(score)
    kd = int(math.ceil(0.1 * n))
    base = float(y.mean())
    out = {"n": n, "n_pos": int(y.sum()), "base_rate": base, "auroc": auroc(score, y)}
    top = y[np.argsort(pos)[:kd]]
    out["top_decile_n"] = kd
    out["top_decile_hit_rate"] = float(top.mean())
    out["top_decile_lift"] = out["top_decile_hit_rate"] / base if base > 0 else np.nan
    for k in KS:
        out[f"precision_at_{k}"] = float(y[np.argsort(pos)[:k]].mean()) if n >= k else np.nan
    out["spearman_count"] = float(_spearman_rows(score[None, :], count[None, :])[0])
    return out


def boot_metrics(score, y, count, geo_ids, idx: np.ndarray) -> dict:
    """Metric distributions over resample index matrix idx (B x n)."""
    s, ys, cs = score[idx], y[idx], count[idx]
    pos = _order_positions(score, geo_ids)[idx]
    n = idx.shape[1]
    out = {"auroc": auroc_from_ranks(rankdata(s, axis=1), ys)}
    base = ys.mean(axis=1)
    kd = int(math.ceil(0.1 * n))
    hit = _topk_rate(pos, ys, kd)
    out["top_decile_hit_rate"] = hit
    with np.errstate(invalid="ignore", divide="ignore"):
        out["top_decile_lift"] = np.where(base > 0, hit / base, np.nan)
    for k in KS:
        out[f"precision_at_{k}"] = _topk_rate(pos, ys, k) if n >= k else np.full(len(idx), np.nan)
    out["spearman_count"] = _spearman_rows(s, cs)
    return out


def boot_index(n: int, B: int, seed: int) -> np.ndarray:
    return np.random.default_rng(seed).integers(0, n, size=(B, n))


def cluster_boot_index(groups: np.ndarray, B: int, seed: int) -> list[np.ndarray]:
    rng = np.random.default_rng(seed)
    ug = np.unique(groups)
    members = {g: np.flatnonzero(groups == g) for g in ug}
    out = []
    for _ in range(B):
        pick = rng.choice(ug, size=len(ug), replace=True)
        out.append(np.concatenate([members[g] for g in pick]))
    return out


def perm_auroc(score: np.ndarray, y: np.ndarray, groups: np.ndarray | None, P: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    r = rankdata(score)
    n = len(y)
    if groups is None:
        Y = np.array([rng.permutation(y) for _ in range(P)])
    else:
        Y = np.empty((P, n), dtype=y.dtype)
        gcodes = pd.factorize(groups)[0]
        for p in range(P):
            key = rng.random(n)
            order = np.lexsort((key, gcodes))          # rows grouped by state, shuffled within
            base = np.lexsort((np.arange(n), gcodes))  # rows grouped by state, original order
            yp = np.empty(n, dtype=y.dtype)
            yp[base] = y[order]
            Y[p] = yp
    return auroc_from_ranks(np.broadcast_to(r, (P, n)), Y)


def pop_stratified_auroc(score, y, pop50) -> float:
    q = pd.qcut(pd.Series(pop50).rank(method="first"), 10, labels=False).to_numpy()
    num = den = 0.0
    for d in range(10):
        m = q == d
        n1 = y[m].sum()
        n0 = m.sum() - n1
        if n1 == 0 or n0 == 0:
            continue
        num += auroc(score[m], y[m]) * n1 * n0
        den += n1 * n0
    return num / den if den > 0 else np.nan


# ------------------------------------------------------------------------------------------------------------------
# evaluation of one pre-period ranking
# ------------------------------------------------------------------------------------------------------------------

def score_columns(frame: pd.DataFrame, oc: pd.DataFrame, rng_seed: int) -> dict[str, np.ndarray]:
    """name -> score aligned to frame rows (composites, components, baselines)."""
    out = {}
    for c in sorted(c for c in frame.columns if c.startswith("composite_") and c != "composite"):
        out[c] = frame[c].to_numpy(float)
    for c in COMPONENT_SCORES:
        out[c] = frame[c].to_numpy(float)
    o = oc.set_index("geo_id").reindex(frame["geo_id"])
    out["baseline_population_within_50km"] = o["population_within_50km"].to_numpy(float)
    out["baseline_population_total"] = o["population_total"].to_numpy(float)
    out["baseline_prior_trials_50km"] = o["prior_trials_pool_n"].to_numpy(float)
    out["baseline_random"] = np.random.default_rng(rng_seed).random(len(frame))
    return out


def kind_of(name: str) -> str:
    if name.startswith("composite_"):
        return "composite"
    if name.startswith("baseline_"):
        return "baseline"
    return SCORE_KIND.get(name, "other")


def evaluate_ranking(frame: pd.DataFrame, oc: pd.DataFrame, key: dict, outcomes_: list[str] | None = None,
                     B: int = N_BOOT, P: int = N_PERM, full: bool = True) -> list[dict]:
    """Every score x outcome: point metrics, county-bootstrap CIs, paired Delta AUROC vs baselines, permutation
    p-values, state-cluster CI, population-stratified AUROC (plan section 5)."""
    elig = frame["rank_eligible"].astype(bool).to_numpy()
    fr = frame[elig].reset_index(drop=True)
    o = oc.set_index("geo_id").reindex(fr["geo_id"]).reset_index()
    scores = score_columns(fr, oc, config.SEED + 7)
    geo_ids = fr["geo_id"].astype(str).to_numpy()
    rows = []
    for ok_name in (outcomes_ or list(OUTCOMES)):
        col, desc = OUTCOMES[ok_name]
        cnt_all = o[col].to_numpy(float)
        y_all = (cnt_all > 0).astype(np.int8)
        pop50_all = o["population_within_50km"].to_numpy(float)
        states_all = o["state_fips"].astype(str).to_numpy()
        # common evaluation set for the paired comparisons: rows where the score and the baselines exist
        base_names = ["baseline_population_within_50km", "composite_burden_only", "baseline_prior_trials_50km"]
        for name, sc in scores.items():
            m = ~np.isnan(sc) & ~np.isnan(cnt_all)
            if m.sum() < 30:
                continue
            s, y, c, g = sc[m], y_all[m], cnt_all[m], geo_ids[m]
            rec = {**key, "outcome": ok_name, "outcome_description": desc, "score": name, "score_kind": kind_of(name),
                   "expected_direction": EXPECTED.get(name, "none")}
            rec.update(point_metrics(s, y, c, g))
            if full and 0 < y.sum() < len(y):
                idx = boot_index(len(s), B, config.SEED + 101)
                bm = boot_metrics(s, y, c, g, idx)
                for k, v in bm.items():
                    rec[f"{k}_ci_low"], rec[f"{k}_ci_high"] = _ci(v)
                for bn in base_names:
                    if bn == name or bn not in scores:
                        continue
                    bsc = scores[bn][m]
                    if np.isnan(bsc).any():
                        continue
                    bb = auroc_from_ranks(rankdata(bsc[idx], axis=1), y[idx])
                    d = bm["auroc"] - bb
                    short = bn.replace("baseline_", "").replace("composite_", "")
                    rec[f"delta_auroc_vs_{short}"] = rec["auroc"] - auroc(bsc, y)
                    rec[f"delta_auroc_vs_{short}_ci_low"], rec[f"delta_auroc_vs_{short}_ci_high"] = _ci(d)
                null = perm_auroc(s, y, None, P, config.SEED + 202)
                rec["perm_p_auroc"] = float((1 + np.sum(null >= rec["auroc"])) / (P + 1))
                rec["perm_p_auroc_lower"] = float((1 + np.sum(null <= rec["auroc"])) / (P + 1))
                rec["perm_null_auroc_p95"] = float(np.percentile(null, 95))
                nullw = perm_auroc(s, y, states_all[m], P, config.SEED + 203)
                rec["perm_within_state_p_auroc"] = float((1 + np.sum(nullw >= rec["auroc"])) / (P + 1))
                rec["perm_within_state_p_auroc_lower"] = float((1 + np.sum(nullw <= rec["auroc"])) / (P + 1))
                # precision@25 permutation p (national)
                rng = np.random.default_rng(config.SEED + 204)
                pos = _order_positions(s, g)
                top25 = np.argsort(pos)[:25]
                nullp = np.array([rng.permutation(y)[top25].mean() for _ in range(P)])
                rec["perm_p_precision_at_25"] = float((1 + np.sum(nullp >= rec["precision_at_25"])) / (P + 1))
                # state-cluster bootstrap (AUROC, Delta vs population)
                cidx = cluster_boot_index(states_all[m], B, config.SEED + 303)
                ca, cd = [], []
                p50 = pop50_all[m]
                for ii in cidx:
                    yy = y[ii]
                    if yy.min() == yy.max():
                        continue
                    a1 = auroc(s[ii], yy)
                    ca.append(a1)
                    if not np.isnan(p50).any():
                        cd.append(a1 - auroc(p50[ii], yy))
                rec["auroc_cluster_ci_low"], rec["auroc_cluster_ci_high"] = _ci(np.array(ca))
                if cd:
                    rec["delta_auroc_vs_population_within_50km_cluster_ci_low"], \
                        rec["delta_auroc_vs_population_within_50km_cluster_ci_high"] = _ci(np.array(cd))
                if not np.isnan(p50).any():
                    rec["auroc_population_stratified"] = pop_stratified_auroc(s, y, p50)
                    ps = [pop_stratified_auroc(s[ii], y[ii], p50[ii]) for ii in idx[:200]]
                    rec["auroc_population_stratified_ci_low"], rec["auroc_population_stratified_ci_high"] = _ci(
                        np.array(ps))
                    rec["partial_spearman_count_given_log_pop50"] = F.partial_spearman(s, c, np.log1p(p50))
            rows.append(rec)
    return rows


def random_baseline(n: int, y: np.ndarray, R: int = 1000) -> dict:
    rng = np.random.default_rng(config.SEED + 404)
    aucs, p = [], {k: [] for k in KS}
    for _ in range(R):
        s = rng.random(n)
        aucs.append(auroc(s, y))
        o = np.argsort(-s)
        for k in KS:
            p[k].append(y[o[:k]].mean())
    out = {"random_auroc_mean": float(np.mean(aucs)), "random_auroc_p025": float(np.percentile(aucs, 2.5)),
           "random_auroc_p975": float(np.percentile(aucs, 97.5))}
    for k in KS:
        out[f"random_precision_at_{k}_mean"] = float(np.mean(p[k]))
        out[f"random_precision_at_{k}_p025"] = float(np.percentile(p[k], 2.5))
        out[f"random_precision_at_{k}_p975"] = float(np.percentile(p[k], 97.5))
    return out


def population_adjusted(frame: pd.DataFrame, oc: pd.DataFrame, key: dict, outcome: str = "O1",
                        names=("composite_equal", "composite_evidence_weighted", "composite_study_partner",
                               "composite_access_gap", *COMPONENT_SCORES), B: int = N_BOOT) -> list[dict]:
    """Logistic O ~ log1p(pop50) + log1p(pop) + z(score) [+ log1p(prior trials 50 km)]: OR per SD (county-bootstrap
    CI), likelihood-ratio p (plan section 5.2 c-d)."""
    import statsmodels.api as sm
    elig = frame["rank_eligible"].astype(bool).to_numpy()
    fr = frame[elig].reset_index(drop=True)
    o = oc.set_index("geo_id").reindex(fr["geo_id"]).reset_index()
    col = OUTCOMES[outcome][0]
    y_all = (o[col].to_numpy(float) > 0).astype(float)
    X0_all = np.column_stack([np.log1p(o["population_within_50km"].to_numpy(float)),
                              np.log1p(o["population_total"].to_numpy(float))])
    prior_all = np.log1p(o["prior_trials_pool_n"].to_numpy(float))
    rows = []
    rng = np.random.default_rng(config.SEED + 505)
    for name in names:
        if name not in fr:
            continue
        sc = fr[name].to_numpy(float)
        m = ~np.isnan(sc) & ~np.isnan(X0_all).any(axis=1)
        y, X0, prior = y_all[m], X0_all[m], prior_all[m]
        z = (sc[m] - sc[m].mean()) / sc[m].std()
        for model, extra in (("population", None), ("population_plus_prior_activity", prior)):
            base = X0 if extra is None else np.column_stack([X0, extra])
            Xr = sm.add_constant(base)
            Xf = sm.add_constant(np.column_stack([base, z]))
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    fr_ = sm.Logit(y, Xf).fit(disp=0, maxiter=200)
                    rr_ = sm.Logit(y, Xr).fit(disp=0, maxiter=200)
            except Exception as e:  # separation etc.: reported, never silently dropped
                rows.append({**key, "outcome": outcome, "score": name, "model": model, "status": f"fit failed: {e}"})
                continue
            lr = 2 * (fr_.llf - rr_.llf)
            ors = []
            for _ in range(B):
                ii = rng.integers(0, len(y), len(y))
                try:
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        b = sm.Logit(y[ii], Xf[ii]).fit(disp=0, maxiter=100).params[-1]
                    ors.append(np.exp(b))
                except Exception:
                    ors.append(np.nan)
            lo, hi = _ci(np.array(ors))
            auc_r = auroc(rr_.predict(Xr), y)
            auc_f = auroc(fr_.predict(Xf), y)
            rows.append({**key, "outcome": outcome, "score": name, "score_kind": kind_of(name),
                         "expected_direction": EXPECTED.get(name, "none"), "model": model, "n": int(len(y)),
                         "n_pos": int(y.sum()), "odds_ratio_per_sd": float(np.exp(fr_.params[-1])),
                         "odds_ratio_ci_low": lo, "odds_ratio_ci_high": hi,
                         "odds_ratio_wald_ci_low": float(np.exp(fr_.conf_int()[-1][0])),
                         "odds_ratio_wald_ci_high": float(np.exp(fr_.conf_int()[-1][1])),
                         "lr_stat": float(lr), "lr_p": float(stats.chi2.sf(lr, 1)),
                         "auroc_model_without_score": auc_r, "auroc_model_with_score": auc_f,
                         "delta_auroc_in_sample": auc_f - auc_r, "status": "ok"})
    return rows


# ------------------------------------------------------------------------------------------------------------------
# state level, burden persistence, stability
# ------------------------------------------------------------------------------------------------------------------

def state_level(frame: pd.DataFrame, os_: pd.DataFrame, key: dict, B: int = N_BOOT) -> list[dict]:
    fr = frame.reset_index(drop=True)
    o = os_.set_index("geo_id").reindex(fr["geo_id"]).reset_index()
    rows = []
    scores = {c: fr[c].to_numpy(float) for c in fr.columns if c.startswith("composite_") and c != "composite"}
    scores.update({c: fr[c].to_numpy(float) for c in COMPONENT_SCORES})
    scores["baseline_population_total"] = o["population_total"].to_numpy(float)
    pop = o["population_total"].to_numpy(float)
    rng = np.random.default_rng(config.SEED + 606)
    idx = rng.integers(0, len(fr), size=(B, len(fr)))
    for oname, col in (("O1_state", "o1_trials_in_state_n"), ("O3_state", "o3_nih_new_core_in_state_n"),
                       ("O4_state", "o4_new_facilities_in_state_n")):
        cnt = o[col].to_numpy(float)
        for per in ("count", "per_100k"):
            v = cnt if per == "count" else cnt / pop * 1e5
            for name, sc in scores.items():
                m = ~np.isnan(sc) & ~np.isnan(v)
                if m.sum() < 10:
                    continue
                rho = float(_spearman_rows(sc[m][None, :], v[m][None, :])[0])
                ii = idx[:, :m.sum()] % m.sum()
                bs = _spearman_rows(sc[m][ii], v[m][ii])
                lo, hi = _ci(bs)
                rows.append({**key, "outcome": oname, "scale": per, "score": name, "score_kind": kind_of(name),
                             "expected_direction": EXPECTED.get(name, "none"), "n_states": int(m.sum()),
                             "n_states_positive": int((cnt[m] > 0).sum()), "spearman": rho, "spearman_ci_low": lo,
                             "spearman_ci_high": hi})
    return rows


def burden_persistence(oc: pd.DataFrame, frame: pd.DataFrame, T: int) -> dict:
    fr = frame[frame["rank_eligible"].astype(bool)]
    o = oc.set_index("geo_id").reindex(fr["geo_id"])
    a, b = o["o5_cms51_T"].to_numpy(float), o["o5_cms51_latest_final"].to_numpy(float)
    m = ~np.isnan(a) & ~np.isnan(b)
    rho = F.spearman(a[m], b[m])
    lo, hi = F.bootstrap_spearman_ci(a[m], b[m])
    d = b[m] - a[m]
    top = pd.Series(a[m]).rank(ascending=False, method="first").to_numpy() <= math.ceil(0.1 * m.sum())
    top_b = pd.Series(b[m]).rank(ascending=False, method="first").to_numpy() <= math.ceil(0.1 * m.sum())
    return {"T": T, "n": int(m.sum()), "spearman_cms51_T_vs_latest_final": rho, "ci_low": lo, "ci_high": hi,
            "mean_change_pp": float(d.mean()), "median_change_pp": float(np.median(d)),
            "spearman_change_vs_T_value": F.spearman(a[m], d),
            "top_decile_retained_share": float((top & top_b).sum() / top.sum()),
            "latest_final_year": CMS_LATEST_FINAL}


def loso(view: View, key: dict) -> list[dict]:
    """Leave-one-state-out normalisation (plan section 7) for the county ranking of the primary query."""
    full = score_view(view)
    geo = view.S.levels["county"].geo
    st = geo["state_fips"].astype(str).to_numpy()
    fc = full["composite_equal"].to_numpy(float)
    elig = full["rank_eligible"].astype(bool).to_numpy()
    full_rank = full["rank_equal"].to_numpy(float)
    rows = []
    for s in sorted(np.unique(st)):
        ref = st != s
        lo = score_view(view, ref=ref)
        lc = lo["composite_equal"].to_numpy(float)
        ins = (st == s) & elig & ~np.isnan(fc) & ~np.isnan(lc)
        oth = (st != s) & elig & ~np.isnan(fc) & ~np.isnan(lc)
        top_f = set(np.flatnonzero(full_rank <= 25))
        top_l = set(np.flatnonzero(lo["rank_equal"].to_numpy(float) <= 25))
        rows.append({**key, "state_fips": s, "state_abbr": str(geo.loc[st == s, "state_abbr"].iloc[0]),
                     "n_state_counties_ranked": int(ins.sum()),
                     "spearman_state_counties": F.spearman(fc[ins], lc[ins]) if ins.sum() >= 3 else np.nan,
                     "max_abs_composite_change_state": float(np.nanmax(np.abs(fc[ins] - lc[ins]))) if ins.any() else np.nan,
                     "spearman_other_counties": F.spearman(fc[oth], lc[oth]),
                     "top25_overlap_national": len(top_f & top_l),
                     "state_counties_in_top25_full": int(sum(1 for i in top_f if st[i] == s)),
                     "state_counties_in_top25_loso": int(sum(1 for i in top_l if st[i] == s))})
    return rows


def rank_agreement(a: pd.DataFrame, b: pd.DataFrame, la: str, lb: str) -> dict:
    x = a.set_index("geo_id")
    y = b.set_index("geo_id")
    common = x.index[x["rank_eligible"].astype(bool)].intersection(y.index[y["rank_eligible"].astype(bool)])
    ca, cb = x.loc[common, "composite_equal"].to_numpy(float), y.loc[common, "composite_equal"].to_numpy(float)
    m = ~np.isnan(ca) & ~np.isnan(cb)
    ra = x.loc[common, "rank_equal"]
    rb = y.loc[common, "rank_equal"]
    out = {"ranking_a": la, "ranking_b": lb, "n_common": int(m.sum()), "spearman": F.spearman(ca[m], cb[m]),
           "kendall_tau_b": float(stats.kendalltau(ca[m], cb[m]).statistic)}
    for n in (10, 25):
        ta, tb = set(x.index[x["rank_equal"] <= n]), set(y.index[y["rank_equal"] <= n])
        out[f"top{n}_overlap"] = len(ta & tb)
    return out


# ------------------------------------------------------------------------------------------------------------------
# run
# ------------------------------------------------------------------------------------------------------------------

KEEP_COLS = ["object_id", "condition_id", "measurement_id", "geo_level", "geo_id", "geo_name", "state_fips",
             "state_abbr", "population_total", "population_within_50km", "rank_eligible", "burden_value",
             "burden_evidence_level", "burden_members_contributing", "burden_incomplete", "burden_pct",
             "svi_overall", "vulnerability_pct", "diagnostic_desert", "diagnostic_desert_pct",
             "cc_impl_providers_in_geo_n", "cc_impl_facilities_in_geo_n", "clinic_capacity_pct",
             "rr_condition_trials_pool_n", "rr_condition_nih_core_pool_n", "rr_tech_trials_pool_n",
             "research_readiness_pct", "reach", "expected_yield", "measurement_evidence", "evidence_weighted_status"]


def _rank_frame(df: pd.DataFrame, T, variant: str, query: str) -> pd.DataFrame:
    cols = [c for c in KEEP_COLS if c in df.columns]
    cols += [c for c in df.columns if (c.startswith("composite_") or c.startswith("rank_")) and c not in cols]
    out = df[cols].copy()
    tl = "current" if T is None else str(T)
    out.insert(0, "holdout_cutoff", tl)
    out.insert(1, "holdout_variant", variant)
    out.insert(2, "holdout_query", query)
    out["object_id"] = [f"temporal_holdout:{tl}|{variant}|{query}|{g}" for g in out["geo_id"]]
    return out


def view_check(n_tol: float = 1e-9) -> pd.DataFrame:
    """T=None view reproduces the stored primary-query composites (county and state) and the stored desert."""
    stored = read_table(O.TABLE, columns=["condition_id", "measurement_id", "geo_level", "geo_id",
                                          "composite_equal", "composite_evidence_weighted", "rank_equal"])
    v = make_view(None)
    rows = []
    for level in O.LEVELS:
        cur = score_view(v, level=level)
        s = stored[(stored["condition_id"] == QUERY) & (stored["measurement_id"] == MEASUREMENT)
                   & (stored["geo_level"] == level)].set_index("geo_id").reindex(cur["geo_id"])
        for c in ("composite_equal", "composite_evidence_weighted", "rank_equal"):
            a, b = cur[c].to_numpy(float), s[c].to_numpy(float)
            both = ~np.isnan(a) & ~np.isnan(b)
            rows.append({"check": f"{level}: view(None) {c} vs stored deployment_opportunities",
                         "n": int(both.sum()), "nan_pattern_equal": bool((np.isnan(a) == np.isnan(b)).all()),
                         "max_abs_diff": float(np.max(np.abs(a[both] - b[both]))) if both.any() else np.nan})
    f0 = O._features()
    for (cid, level), raw in v.raw.items():
        st = f0[(f0["condition_id"] == cid) & (f0["geo_level"] == level)].set_index("geo_id").reindex(raw["geo_id"])
        mine = v.features[(v.features["condition_id"] == cid) & (v.features["geo_level"] == level)]
        a, b = mine["diagnostic_desert"].to_numpy(float), st["diagnostic_desert"].to_numpy(float)
        both = ~np.isnan(a) & ~np.isnan(b)
        rows.append({"check": f"{level} {cid}: rebuilt desert vs stored geo_condition_features", "n": int(both.sum()),
                     "nan_pattern_equal": bool((np.isnan(a) == np.isnan(b)).all()),
                     "max_abs_diff": float(np.max(np.abs(a[both] - b[both]))) if both.any() else np.nan})
        with engine_view(v):
            _, rate = O.relevant_provider_counts(cid, level)
        a, b = rate, st["relevant_providers_per_100k"].to_numpy(float)
        both = ~np.isnan(a) & ~np.isnan(b)
        rows.append({"check": f"{level} {cid}: engine relevant_provider_counts vs stored relevant_providers_per_100k",
                     "n": int(both.sum()), "nan_pattern_equal": bool((np.isnan(a) == np.isnan(b)).all()),
                     "max_abs_diff": float(np.max(np.abs(a[both] - b[both]))) if both.any() else np.nan})
    return pd.DataFrame(rows)


def _w(df: pd.DataFrame, name: str) -> str:
    TABLES.mkdir(parents=True, exist_ok=True)
    p = TABLES / f"temporal_holdout_{name}.csv"
    df.to_csv(p, index=False)
    return str(p.relative_to(config.PROJECT_ROOT))


def run(n_boot: int = N_BOOT, n_perm: int = N_PERM, cutoffs=CUTOFFS, write: bool = True) -> dict:
    t0 = time.time()
    built = utc_now_iso()
    checks = view_check()
    log(f"view check done ({time.time() - t0:.0f}s)")
    rank_frames, metrics, popadj, state_rows, leak, burden_rows, rand_rows, top_rows = [], [], [], [], [], [], [], []
    info_rows, loso_rows, rolling = [], [], []
    current = {lvl: score_view(make_view(None), level=lvl) for lvl in O.LEVELS}
    rank_frames += [_rank_frame(current[l], None, "current", QUERY) for l in O.LEVELS]
    primary_frames = {}
    for T in cutoffs:
        oc, os_, info = outcomes(T)
        info_rows.append({"T": T, **info})
        for variant in VARIANTS:
            if variant != "A_clean" and T != PRIMARY_T:
                continue
            view = make_view(T, variant)
            leak.append({"T": T, "variant": variant, **view.stats})
            queries = [(QUERY, "county"), (QUERY, "state")]
            if variant == "A_clean":
                queries.append(("me_cfs", "county"))
            for q, level in queries:
                fr = score_view(view, condition=q, level=level)
                rank_frames.append(_rank_frame(fr, T, variant, q))
                key = {"T": T, "variant": variant, "clean_holdout": VARIANTS[variant]["clean"], "query": q,
                       "level": level}
                if level == "state":
                    state_rows += state_level(fr, os_, key, B=n_boot)
                    continue
                if variant == "A_clean" and q == QUERY:
                    primary_frames[T] = fr
                full = (T == PRIMARY_T) or variant == "A_clean"
                outs = list(OUTCOMES) if (variant == "A_clean" and q == QUERY) else ["O1", "O2", "O4"]
                metrics += evaluate_ranking(fr, oc, key, outs, B=n_boot, P=n_perm, full=full)
                log(f"T={T} {variant} {q} {level}: evaluated ({time.time() - t0:.0f}s)")
                if T == PRIMARY_T and variant == "A_clean" and q == QUERY:
                    popadj += population_adjusted(fr, oc, key, "O1", B=n_boot)
                    popadj += population_adjusted(fr, oc, key, "O4", names=("composite_equal", "research_readiness_pct",
                                                                           "diagnostic_desert_pct"), B=n_boot)
                    elig = fr[fr["rank_eligible"].astype(bool)]
                    o = oc.set_index("geo_id").reindex(elig["geo_id"])
                    for oname in ("O1", "O4", "O3"):
                        y = (o[OUTCOMES[oname][0]].to_numpy(float) > 0).astype(float)
                        rand_rows.append({**key, "outcome": oname, "n": len(y), "base_rate": float(y.mean()),
                                          **random_baseline(len(y), y)})
                    top = elig.sort_values("rank_equal").head(25)
                    oo = oc.set_index("geo_id").reindex(top["geo_id"])
                    top_rows.append(pd.DataFrame({
                        "T": T, "rank_equal_pre_T": top["rank_equal"].to_numpy(), "geo_id": top["geo_id"].to_numpy(),
                        "geo_name": top["geo_name"].to_numpy(),
                        "research_readiness_pct": top["research_readiness_pct"].to_numpy(),
                        "diagnostic_desert_pct": top["diagnostic_desert_pct"].to_numpy(),
                        "prior_trials_pool_n": oo["prior_trials_pool_n"].to_numpy(),
                        "post_trials_pool_n": oo["o1_trials_pool_n"].to_numpy(),
                        "post_new_facilities_pool_n": oo["o4_new_facilities_pool_n"].to_numpy(),
                        "post_nih_new_core_in_county_n": oo["o3_nih_new_core_in_county_n"].to_numpy(),
                        "current_rank_equal": current["county"].set_index("geo_id")["rank_equal"].reindex(
                            top["geo_id"]).to_numpy()}))
                    burden_rows.append(burden_persistence(oc, fr, T))
                    _DECILES["table"] = compute_deciles(fr, oc)
                    # outcomes attached to the processed ranking rows
                    rank_frames[-1] = rank_frames[-1].merge(
                        oc[["geo_id", "prior_trials_pool_n", "o1_trials_pool_n", "o1_trials_in_county_n",
                            "o2_trials_pool_n", "o3_nih_new_core_in_county_n", "o4_new_facilities_pool_n"]].rename(
                            columns=lambda c: c if c == "geo_id" else f"post_{c}" if not c.startswith("prior") else c),
                        on="geo_id", how="left")
        if T != PRIMARY_T:
            burden_rows.append(burden_persistence(oc, primary_frames[T], T))
    # stability
    log(f"LOSO ({time.time() - t0:.0f}s)")
    loso_rows += loso(make_view(None), {"ranking": "current"})
    loso_rows += loso(make_view(PRIMARY_T, "A_clean"), {"ranking": f"T={PRIMARY_T} A_clean"})
    labels = {f"T={T}": primary_frames[T] for T in sorted(primary_frames)}
    labels["current"] = current["county"]
    names = list(labels)
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            rolling.append(rank_agreement(labels[a], labels[b], a, b))
    out = {"metrics": pd.DataFrame(metrics), "population_adjusted": pd.DataFrame(popadj),
           "state": pd.DataFrame(state_rows), "leakage": pd.DataFrame(leak), "burden_persistence": pd.DataFrame(burden_rows),
           "random_baseline": pd.DataFrame(rand_rows), "top25": pd.concat(top_rows, ignore_index=True),
           "outcome_counts": pd.DataFrame(info_rows), "loso": pd.DataFrame(loso_rows), "rolling": pd.DataFrame(rolling),
           "view_check": checks}
    out["loso_summary"] = loso_summary(out["loso"])
    out["deciles"] = _DECILES.get("table", pd.DataFrame())
    if write:
        paths = {k: _w(v, k) for k, v in out.items()}
        rk = pd.concat(rank_frames, ignore_index=True)
        rk = add_provenance(
            rk, data_layer="derived",
            source_name=f"{PRODUCER} (opportunity engine on dated input views: clinical_trials, trial_sites, "
                        "trial_conditions, facility_trials, nih_projects, facility_nih_projects, providers, facilities, "
                        "facilities__hrsa, geo_condition_burden (CMS MMD 51), geo_condition_features)",
            source_version=f"temporal holdout {built}; plan {PLAN}", retrieved_at=built,
            evidence_type="derived_score", source_record_id="object_id",
            source_geographic_resolution="county",
            evidence_level=rk["burden_evidence_level"].astype(str).to_numpy(),
            provenance_notes=("Pre-period deployment rankings rebuilt with inputs dated on or before the cutoff "
                              "(holdout_cutoff; 'current' = today's inputs) by the unchanged scoring engine, with "
                              "post-period research-activity outcomes (post_*). " + GUARDRAIL))
        rk["source_geographic_resolution"] = np.where(rk["geo_level"] == "state", "state", "county")
        assert rk["object_id"].is_unique
        write_table(rk, TABLE, producer=PRODUCER,
                    description="Temporal holdout (Test 8): deployment rankings rebuilt from pre-cutoff inputs, every "
                                "component and composite kept, with post-cutoff research-activity outcomes. "
                                "docs/ANALYSIS_PLAN_TEMPORAL_HOLDOUT.md")
        fig = draft_figure(out)
        write_results_md(out, paths, fig, built)
    log(f"done ({time.time() - t0:.0f}s)")
    return out


def loso_summary(l: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for r, g in l.groupby("ranking"):
        s = g["spearman_state_counties"]
        rows.append({"ranking": r, "n_states": int(len(g)), "n_states_with_3plus_ranked": int(s.notna().sum()),
                     "median_spearman_state_counties": float(s.median()), "min_spearman_state_counties": float(s.min()),
                     "states_below_0_8": SEP.join(g.loc[s < 0.8, "state_abbr"].astype(str)
                                                  if "state_abbr" in g else g.loc[s < 0.8, "state_fips"].astype(str)),
                     "median_spearman_other_counties": float(g["spearman_other_counties"].median()),
                     "min_spearman_other_counties": float(g["spearman_other_counties"].min()),
                     "median_top25_overlap": float(g["top25_overlap_national"].median()),
                     "min_top25_overlap": int(g["top25_overlap_national"].min())})
    return pd.DataFrame(rows)


# ------------------------------------------------------------------------------------------------------------------
# figure and results markdown
# ------------------------------------------------------------------------------------------------------------------

def draft_figure(out: dict) -> str:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    m = out["metrics"]
    p = m[(m["T"] == PRIMARY_T) & (m["variant"] == "A_clean") & (m["query"] == QUERY) & (m["outcome"] == "O1")]
    order = ["composite_equal", "composite_evidence_weighted", "composite_study_partner", "composite_capacity_led",
             "composite_access_gap", "research_readiness_pct", "clinic_capacity_pct", "burden_pct",
             "vulnerability_pct", "diagnostic_desert_pct", "baseline_prior_trials_50km",
             "baseline_population_within_50km", "baseline_population_total", "composite_burden_only",
             "baseline_random"]
    p = p.set_index("score").reindex([o for o in order if o in set(p["score"])])
    fig, ax = plt.subplots(1, 2, figsize=(13, 6))
    colors = {"composite": "#3b6ea5", "component": "#7a9a3a", "baseline": "#8a8a8a"}
    yv = np.arange(len(p))[::-1]
    for i, (name, r) in enumerate(p.iterrows()):
        kind = "baseline" if name == "composite_burden_only" else r["score_kind"]
        ax[0].errorbar(r["auroc"], yv[i], xerr=[[r["auroc"] - r["auroc_ci_low"]], [r["auroc_ci_high"] - r["auroc"]]],
                       fmt="o", color=colors.get(kind, "k"), capsize=3)
    ax[0].axvline(0.5, color="k", lw=0.8, ls="--")
    ax[0].set_yticks(yv)
    ax[0].set_yticklabels([s.replace("composite_", "composite: ").replace("baseline_", "baseline: ") for s in p.index],
                          fontsize=8)
    ax[0].set_xlabel("AUROC (95% county-bootstrap CI)")
    ax[0].set_title(f"Pre-{PRIMARY_T + 1} ranking vs any new Long COVID / ME/CFS trial\nwithin 50 km after "
                    f"{PRIMARY_T} (O1)", fontsize=10)
    # panel B: hit rate by decile of four scores
    dec = out.get("deciles", pd.DataFrame(columns=["score", "decile", "hit_rate", "base_rate"]))
    for name, lab, c in (("composite_equal", "equal composite", "#3b6ea5"),
                         ("research_readiness_pct", "research readiness", "#7a9a3a"),
                         ("diagnostic_desert_pct", "diagnostic desert", "#c0504d"),
                         ("baseline_population_within_50km", "population within 50 km", "#8a8a8a")):
        d = dec[dec["score"] == name]
        ax[1].plot(d["decile"], d["hit_rate"], marker="o", label=lab, color=c)
    if len(dec):
        ax[1].axhline(dec["base_rate"].iloc[0], color="k", lw=0.8, ls="--", label="base rate")
    ax[1].set_xlabel("decile of the pre-period score (10 = highest)")
    ax[1].set_ylabel("share of counties with a new relevant trial within 50 km")
    ax[1].legend(fontsize=8)
    ax[1].set_title("Hit rate by score decile (eligible counties)", fontsize=10)
    fig.suptitle("DRAFT - Test 8 temporal holdout (proxy outcome: research activity, not clinical value)", fontsize=10)
    fig.tight_layout()
    d = FIGURES / "drafts"
    d.mkdir(parents=True, exist_ok=True)
    path = d / "temporal_holdout_auroc_deciles.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return str(path.relative_to(config.PROJECT_ROOT))


_DECILES: dict = {}


def compute_deciles(frame: pd.DataFrame, oc: pd.DataFrame) -> pd.DataFrame:
    fr = frame[frame["rank_eligible"].astype(bool)].reset_index(drop=True)
    o = oc.set_index("geo_id").reindex(fr["geo_id"])
    y = (o["o1_trials_pool_n"].to_numpy(float) > 0).astype(float)
    rows = []
    sc = score_columns(fr, oc, 0)
    for name in ("composite_equal", "research_readiness_pct", "diagnostic_desert_pct", "baseline_population_within_50km",
                 "composite_evidence_weighted", "clinic_capacity_pct", "vulnerability_pct", "burden_pct"):
        s = sc[name]
        m = ~np.isnan(s)
        q = pd.qcut(pd.Series(s[m]).rank(method="first"), 10, labels=False).to_numpy() + 1
        for d in range(1, 11):
            rows.append({"score": name, "decile": d, "n": int((q == d).sum()), "hit_rate": float(y[m][q == d].mean()),
                         "base_rate": float(y[m].mean())})
    return pd.DataFrame(rows)


def _f(x, nd=3) -> str:
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return "NA"
    if isinstance(x, (int, np.integer)):
        return f"{x:,}"
    return f"{x:.{nd}f}"


def _ci_s(r, k, nd=3) -> str:
    return f"{_f(r.get(k), nd)} [{_f(r.get(k + '_ci_low'), nd)}, {_f(r.get(k + '_ci_high'), nd)}]"


def write_results_md(out: dict, paths: dict, fig: str, built: str) -> None:
    meta = {"built": built, "figure": fig, "tables": paths, "plan": PLAN, "n_boot": N_BOOT, "n_perm": N_PERM}
    (TABLES / "temporal_holdout_run_metadata.json").write_text(json.dumps(meta, indent=1))
    report()


SCORE_LABEL = {
    "composite_equal": "composite, equal weights (primary)", "composite_evidence_weighted": "composite, evidence_weighted",
    "composite_study_partner": "composite, study_partner", "composite_capacity_led": "composite, capacity_led",
    "composite_burden_led": "composite, burden_led", "composite_equity_led": "composite, equity_led",
    "composite_access_gap": "composite, access_gap", "composite_saturation_adjusted": "composite, saturation_adjusted",
    "composite_burden_only": "baseline: burden only (engine rank)", "research_readiness_pct": "component: research_readiness",
    "clinic_capacity_pct": "component: clinic_capacity", "burden_pct": "component: burden (ME/CFS CMS proxy)",
    "vulnerability_pct": "component: vulnerability (SVI)", "diagnostic_desert_pct": "component: diagnostic_desert",
    "baseline_population_within_50km": "baseline: population within 50 km",
    "baseline_population_total": "baseline: county population", "baseline_prior_trials_50km":
    "baseline: prior relevant trials within 50 km", "baseline_random": "baseline: one random ranking",
}
PRIMARY_ORDER = ["composite_equal", "composite_evidence_weighted", "composite_study_partner", "composite_capacity_led",
                 "composite_burden_led", "composite_equity_led", "composite_access_gap", "composite_saturation_adjusted",
                 "research_readiness_pct", "clinic_capacity_pct", "burden_pct", "vulnerability_pct",
                 "diagnostic_desert_pct", "baseline_prior_trials_50km", "baseline_population_within_50km",
                 "baseline_population_total", "composite_burden_only", "baseline_random"]


def _read(name: str) -> pd.DataFrame:
    p = TABLES / f"temporal_holdout_{name}.csv"
    return pd.read_csv(p, dtype={"state_fips": str, "geo_id": str}) if p.exists() else pd.DataFrame()


def composite_reading(r) -> str:
    """Pre-specified rule (plan 5.3): (a) AUROC CI lower > 0.5; (b) permutation p < 0.05; (c) Delta vs population
    within 50 km CI lower > 0."""
    a = r["auroc_ci_low"] > 0.5
    b = r["perm_p_auroc"] < 0.05
    c = r.get("delta_auroc_vs_population_within_50km_ci_low", np.nan) > 0
    if a and b and c:
        return "passes (a, b, c)"
    if a and b:
        return "explained by population: (a) and (b) hold, (c) fails"
    if r["auroc_ci_high"] < 0.5:
        return "fails: AUROC below 0.5 (points away from later activity)"
    return "fails: AUROC not above 0.5"


def component_reading(name: str, r) -> str:
    exp = EXPECTED.get(name, "none")
    lo, hi = r["auroc_ci_low"], r["auroc_ci_high"]
    if name == "diagnostic_desert_pct" or name in ("composite_access_gap", "composite_equity_led"):
        gp = r.get("top_decile_lift_ci_high", np.nan) < 1
        return ("gap persistence holds (top-decile lift CI upper < 1)" if gp else
                "gap persistence NOT shown (top-decile lift CI reaches 1)") + \
            (f"; AUROC below 0.5 as expected" if hi < 0.5 else "; AUROC not below 0.5")
    if exp == "positive":
        sup = lo > 0.5 and r.get("delta_auroc_vs_population_within_50km_ci_low", np.nan) > 0
        if name.startswith("baseline"):
            return "as expected" if lo > 0.5 else "not above 0.5"
        return ("supported (AUROC > 0.5 and beats population)" if sup else
                "above 0.5 but does not beat population" if lo > 0.5 else "not above 0.5")
    if exp == "negative":
        return "negative, as expected" if hi < 0.5 else ("positive, against expectation" if lo > 0.5 else
                                                         "indistinguishable from 0.5")
    return ("above 0.5" if lo > 0.5 else "below 0.5" if hi < 0.5 else "indistinguishable from 0.5") + \
        " (no directional prediction)"


def report() -> str:
    """results/TEMPORAL_HOLDOUT_RESULTS.md from the tables written by run() (numbers only from those tables)."""
    import datetime as _dt
    meta_p = TABLES / "temporal_holdout_run_metadata.json"
    meta = json.loads(meta_p.read_text()) if meta_p.exists() else {}
    m = _read("metrics")
    pa = _read("population_adjusted")
    st = _read("state")
    lk = _read("leakage")
    bp = _read("burden_persistence")
    rb = _read("random_baseline")
    top = _read("top25")
    oc = _read("outcome_counts")
    ls = _read("loso_summary")
    lo = _read("loso")
    ro = _read("rolling")
    vc = _read("view_check")
    dec = _read("deciles")
    L = []
    w = L.append
    sel = lambda T, v, q, o: m[(m["T"] == T) & (m["variant"] == v) & (m["query"] == q) & (m["outcome"] == o)] \
        .set_index("score")  # noqa: E731
    P = sel(PRIMARY_T, "A_clean", QUERY, "O1")
    eq, ew = P.loc["composite_equal"], P.loc["composite_evidence_weighted"]
    rr, cc, ds = P.loc["research_readiness_pct"], P.loc["clinic_capacity_pct"], P.loc["diagnostic_desert_pct"]
    pop, prior = P.loc["baseline_population_within_50km"], P.loc["baseline_prior_trials_50km"]
    bo = P.loc["composite_burden_only"]
    best = P.drop(index=[i for i in P.index if i.startswith("baseline")]).sort_values("auroc", ascending=False)
    w("# Temporal holdout (Test 8): would the engine have flagged, before the fact, where Long COVID / ME/CFS "
      "research activity appeared?")
    w("")
    w(f"_Generated by `{PRODUCER}` at {meta.get('built', UNKNOWN)} from the tables it wrote in the same run "
      "(`results/tables/temporal_holdout_*.csv`, processed table `temporal_holdout_rankings`). Plan, pre-specified "
      f"before any holdout outcome was computed: `{PLAN}`. Bootstrap resamples {meta.get('n_boot', N_BOOT):,}, "
      f"permutations {meta.get('n_perm', N_PERM):,}, seed {config.SEED}._")
    w("")
    w("**Read this first.** The outcome is a **proxy outcome**: registered trials, new trial sites and NIH awards that "
      "appeared after the cutoff. A county that later attracted research activity was *research-ready*; that is not "
      "evidence that deploying a measurement there would be clinically valuable, detect more cases or benefit "
      "patients, and nothing here measures that. The desert side of the engine is built to point at places with "
      "little activity, so it is expected NOT to predict activity; for it the checkable statement is gap "
      "persistence. Everything is ecological (places, registries, facilities); nothing describes a person.")
    w("")
    # ---------------------------------------------------------------- headline
    w("## 0. Headline")
    w("")
    w(f"Primary: Long COVID or ME/CFS x wearable autonomic/activity monitoring, county level, cutoff T = {PRIMARY_T} "
      f"(inputs dated <= {PRIMARY_T}-12-31), clean variant A, outcome O1 = any Long COVID or ME/CFS trial first "
      f"posted after T with a U.S. site in the county or within 50 km. {int(eq['n']):,} ranking-eligible counties, "
      f"{int(eq['n_pos']):,} positive (base rate {_f(eq['base_rate'])}).")
    w("")
    w(f"* **Equal-weight composite (the default ranking): {composite_reading(eq)}.** AUROC {_ci_s(eq, 'auroc')} "
      f"(state-cluster CI {_f(eq['auroc_cluster_ci_low'])}-{_f(eq['auroc_cluster_ci_high'])}; permutation p "
      f"{_f(eq['perm_p_auroc'])}); population within 50 km alone reaches {_ci_s(pop, 'auroc')}, so Delta AUROC vs "
      f"population = {_ci_s(eq, 'delta_auroc_vs_population_within_50km')}. Top-decile lift {_ci_s(eq, 'top_decile_lift', 2)}; "
      f"precision@25 {_f(eq['precision_at_25'], 2)} vs base rate {_f(eq['base_rate'], 2)}.")
    w(f"* **evidence_weighted composite: {composite_reading(ew)}.** AUROC {_ci_s(ew, 'auroc')}; Delta vs population "
      f"{_ci_s(ew, 'delta_auroc_vs_population_within_50km')}.")
    w(f"* **Research readiness (the component that should predict activity): {component_reading('research_readiness_pct', rr)}.** "
      f"AUROC {_ci_s(rr, 'auroc')}; Delta vs population {_ci_s(rr, 'delta_auroc_vs_population_within_50km')}; "
      f"Delta vs prior activity alone {_ci_s(rr, 'delta_auroc_vs_prior_trials_50km')}; population-stratified AUROC "
      f"{_ci_s(rr, 'auroc_population_stratified')}.")
    w(f"* **Diagnostic desert (built to point where activity is absent): {component_reading('diagnostic_desert_pct', ds)}.** "
      f"AUROC {_ci_s(ds, 'auroc')}; top-decile hit rate {_f(ds['top_decile_hit_rate'], 3)} vs base "
      f"{_f(ds['base_rate'], 3)} (lift {_ci_s(ds, 'top_decile_lift', 2)}).")
    w(f"* **Population is the strongest single predictor** (AUROC {_f(pop['auroc'])}); prior activity alone "
      f"{_f(prior['auroc'])}; the best engine score is {SCORE_LABEL.get(best.index[0], best.index[0])} "
      f"({_f(best['auroc'].iloc[0])}). Burden-only {_f(bo['auroc'])}.")
    padj = pa[(pa["model"] == "population") & (pa["outcome"] == "O1")].set_index("score") if len(pa) else pd.DataFrame()
    if len(padj) and "composite_equal" in padj.index:
        r = padj.loc["composite_equal"]
        r2 = pa[(pa["model"] == "population_plus_prior_activity") & (pa["outcome"] == "O1")].set_index("score")
        w(f"* **Population-adjusted (logistic, log population + log population within 50 km):** equal composite odds "
          f"ratio per SD {_f(r['odds_ratio_per_sd'], 2)} [{_f(r['odds_ratio_ci_low'], 2)}, "
          f"{_f(r['odds_ratio_ci_high'], 2)}] (LR p {_f(r['lr_p'], 4)}); research readiness "
          f"{_f(padj.loc['research_readiness_pct', 'odds_ratio_per_sd'], 2)} "
          f"[{_f(padj.loc['research_readiness_pct', 'odds_ratio_ci_low'], 2)}, "
          f"{_f(padj.loc['research_readiness_pct', 'odds_ratio_ci_high'], 2)}]; adding prior activity to the model, "
          f"equal composite {_f(r2.loc['composite_equal', 'odds_ratio_per_sd'], 2)} "
          f"[{_f(r2.loc['composite_equal', 'odds_ratio_ci_low'], 2)}, {_f(r2.loc['composite_equal', 'odds_ratio_ci_high'], 2)}].")
    if len(ls):
        for _, r in ls.iterrows():
            w(f"* **Leave-one-state-out normalisation ({r['ranking']}):** median within-state Spearman "
              f"{_f(r['median_spearman_state_counties'])} (min {_f(r['min_spearman_state_counties'])}; states < 0.8: "
              f"{r['states_below_0_8'] if isinstance(r['states_below_0_8'], str) and r['states_below_0_8'] else 'none'}); "
              f"other counties median {_f(r['median_spearman_other_counties'], 4)}; national top-25 overlap median "
              f"{_f(r['median_top25_overlap'], 1)} (min {int(r['min_top25_overlap'])}).")
    if len(ro):
        r = ro[(ro["ranking_a"] == "T=2019") & (ro["ranking_b"] == f"T={PRIMARY_T}")]
        rc = ro[(ro["ranking_a"] == f"T={PRIMARY_T}") & (ro["ranking_b"] == "current")]
        if len(r) and len(rc):
            r, rc = r.iloc[0], rc.iloc[0]
            w(f"* **Rolling check:** T=2019 vs T={PRIMARY_T} Spearman {_f(r['spearman'])}, top-25 overlap "
              f"{int(r['top25_overlap'])}; T={PRIMARY_T} vs current Spearman {_f(rc['spearman'])}, top-25 overlap "
              f"{int(rc['top25_overlap'])}, top-10 {int(rc['top10_overlap'])}.")
    w("")
    engine = [i for i in P.index if not i.startswith("baseline") and i != "composite_burden_only"]
    beat_pop = [i for i in engine if P.loc[i, "delta_auroc_vs_population_within_50km_ci_low"] > 0]
    beat_prior = [i for i in engine if P.loc[i, "delta_auroc_vs_prior_trials_50km_ci_low"] > 0]
    adj_pos = [i for i in padj.index if padj.loc[i, "odds_ratio_ci_low"] > 1] if len(padj) else []
    adj_neg = [i for i in padj.index if padj.loc[i, "odds_ratio_ci_high"] < 1] if len(padj) else []
    r2 = pa[(pa["model"] == "population_plus_prior_activity") & (pa["outcome"] == "O1")].set_index("score") \
        if len(pa) else pd.DataFrame()
    adj2_pos = [i for i in r2.index if r2.loc[i, "odds_ratio_ci_low"] > 1] if len(r2) else []
    lab = lambda xs: ", ".join(SCORE_LABEL.get(x, x) for x in xs) or "none"  # noqa: E731
    w("**Pre-specified readings, in one place.**")
    w("")
    w(f"* **Not validated (null against the population baseline):** no engine score beats population within 50 km on "
      f"AUROC (scores whose Delta AUROC CI is above 0: {lab(beat_pop)}). By the pre-specified rule (plan 5.3 (c)) the "
      f"equal-weight composite, and every other composite, is **explained by population**: big, dense places attract "
      f"trials, and a ranking that knows only population does better than any ranking the engine builds. The "
      f"evidence_weighted composite does not rank later activity above chance at all.")
    w(f"* **Partly validated (secondary population control, plan 5.2 (c)-(d)):** holding population fixed, the scores "
      f"whose odds ratio per SD is above 1 are: {lab(adj_pos)}; with prior activity also in the model: "
      f"{lab(adj2_pos)}. The equal-weight composite is not among them. Research readiness also beats the "
      f"prior-activity baseline on AUROC (scores that do: {lab(beat_prior)}). So the readiness side of the engine "
      f"carries information about where research then appeared beyond population and beyond simple persistence, "
      f"but the default composite dilutes it with components that point the other way.")
    w(f"* **Gap persistence holds (the desert side):** the diagnostic desert's top decile had a later-activity rate "
      f"of {_f(ds['top_decile_hit_rate'], 3)} vs {_f(ds['base_rate'], 3)} overall; after the population adjustment the "
      f"scores with odds ratio below 1 are {lab(adj_neg)}. Flagged deserts stayed without new activity; this does not "
      f"show that a deployment there would be valuable.")
    bpct = P.loc["burden_pct"]
    w(f"* **Burden:** the ME/CFS claims proxy (no directional prediction) has AUROC {_ci_s(bpct, 'auroc')}: "
      + ("unrelated to later research activity." if bpct["auroc_ci_low"] <= 0.5 <= bpct["auroc_ci_high"] else
         "related to later research activity in the direction shown."))
    w("* **Cannot be validated here:** clinical value, diagnostic yield, patient benefit, or whether the flagged "
      "places are the right ones for a deployment (section 11).")
    w("")
    # ---------------------------------------------------------------- method
    w("## 1. What was done")
    w("")
    w("* **Same engine, dated inputs.** `opportunity.make_combo` / `combo_frame` and `metric_link.reach_for` ran "
      "unchanged on a view in which every input dated after the cutoff was removed (trials by first-posted date, NIH "
      "rows by fiscal year, NPIs by enumeration date, facilities by their earliest dated source, CMS MMD condition 51 "
      "of claims year T). The same view with no cutoff reproduces the stored `deployment_opportunities` composites "
      "exactly:")
    w("")
    if len(vc):
        w("| check | n | NaN pattern equal | max abs diff |")
        w("|---|---|---|---|")
        for _, r in vc.iterrows():
            w(f"| {r['check']} | {int(r['n']):,} | {r['nan_pattern_equal']} | {_f(r['max_abs_diff'], 12)} |")
    w("")
    w("* **Variants** (plan section 3): A_clean = Long COVID has no burden measure before mid-2022, so its burden is "
      "undefined and the set burden is the ME/CFS CMS proxy of year T (the engine's rule for a member without a "
      "defined burden); B_lc_post_burden = Long COVID burden held at its current HPS value (**NOT a clean holdout**); "
      "A_static_context = A with providers and facilities not date-filtered (leakage bound).")
    w("* **Cutoffs:** T = 2021 (primary), 2020 and 2019 (sensitivity; 2019 is pre-COVID, so all Long COVID activity "
      "is post-period).")
    w("")
    # ---------------------------------------------------------------- outcomes
    w("## 2. Post-period outcomes (counts, all U.S.)")
    w("")
    if len(oc):
        w("| T | relevant trials first posted <= T | relevant trials first posted > T (O1) | ME/CFS trials > T (O2) | "
          "relevant trials starting > T (O1_start) | new relevant NIH core projects (O3) | new facilities (O4) | "
          "facilities with a pre-T relevant trial |")
        w("|---|---|---|---|---|---|---|---|")
        for _, r in oc.iterrows():
            w(f"| {int(r['T'])} | {int(r['prior_relevant_trials']):,} | {int(r['o1_relevant_trials']):,} | "
              f"{int(r['o2_relevant_trials']):,} | {int(r['o1start_relevant_trials']):,} | "
              f"{int(r['o3_nih_new_core_projects']):,} | {int(r['new_facilities']):,} | "
              f"{int(r['facilities_with_pre_T_relevant_trial']):,} |")
    w("")
    w("Base rates among ranking-eligible counties (primary ranking of each T):")
    w("")
    w("| T | outcome | description | counties | positive | base rate |")
    w("|---|---|---|---|---|---|")
    for T in sorted(m["T"].unique(), reverse=True):
        for o in OUTCOMES:
            s_ = sel(T, "A_clean", QUERY, o)
            if "composite_equal" in s_.index:
                r = s_.loc["composite_equal"]
                w(f"| {T} | {o} | {OUTCOMES[o][1]} | {int(r['n']):,} | {int(r['n_pos']):,} | {_f(r['base_rate'])} |")
    w("")
    # ---------------------------------------------------------------- primary table
    w(f"## 3. Primary result: T = {PRIMARY_T}, clean variant, county, O1 (every score, positive and null alike)")
    w("")
    w("AUROC: higher score = predicted to see a new relevant trial within 50 km. CI = county bootstrap; cluster CI = "
      "state-cluster bootstrap; perm p = one-sided (AUROC >= observed) national permutation; within-state p = "
      "permutation within state; pop-strat = AUROC within deciles of population within 50 km.")
    w("")
    w("| score | expected | AUROC [95% CI] | cluster CI | perm p / within-state p | Delta vs pop 50 km | Delta vs "
      "prior activity | pop-strat AUROC | top-decile lift | P@10 | P@25 | P@100 | Spearman (count) | pre-specified reading |")
    w("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for name in [s for s in PRIMARY_ORDER if s in P.index]:
        r = P.loc[name]
        if name.startswith("composite") and name not in ("composite_burden_only", "composite_access_gap",
                                                         "composite_equity_led"):
            reading = composite_reading(r)
        else:
            reading = component_reading(name, r)
        w(f"| {SCORE_LABEL.get(name, name)} | {EXPECTED.get(name, 'none')} | {_ci_s(r, 'auroc')} | "
          f"{_f(r.get('auroc_cluster_ci_low'))}-{_f(r.get('auroc_cluster_ci_high'))} | "
          f"{_f(r.get('perm_p_auroc'))} / {_f(r.get('perm_within_state_p_auroc'))} | "
          f"{_ci_s(r, 'delta_auroc_vs_population_within_50km') if name != 'baseline_population_within_50km' else '-'} | "
          f"{_ci_s(r, 'delta_auroc_vs_prior_trials_50km') if name != 'baseline_prior_trials_50km' else '-'} | "
          f"{_f(r.get('auroc_population_stratified'))} | {_ci_s(r, 'top_decile_lift', 2)} | "
          f"{_f(r['precision_at_10'], 2)} | {_f(r['precision_at_25'], 2)} | {_f(r['precision_at_100'], 2)} | "
          f"{_ci_s(r, 'spearman_count')} | {reading} |")
    w("")
    if len(rb):
        r = rb[(rb["T"] == PRIMARY_T) & (rb["outcome"] == "O1")]
        if len(r):
            r = r.iloc[0]
            w(f"Random baseline (1,000 random rankings of the same counties): AUROC mean {_f(r['random_auroc_mean'])} "
              f"(95% range {_f(r['random_auroc_p025'])}-{_f(r['random_auroc_p975'])}); precision@25 mean "
              f"{_f(r['random_precision_at_25_mean'], 3)} (95% range {_f(r['random_precision_at_25_p025'], 2)}-"
              f"{_f(r['random_precision_at_25_p975'], 2)}); precision@100 mean {_f(r['random_precision_at_100_mean'], 3)}. "
              "The single-draw 'one random ranking' row above is one such draw.")
            w("")
    if len(dec):
        w("Hit rate by decile of the pre-period score (1 = lowest, 10 = highest; `temporal_holdout_deciles.csv`):")
        w("")
        piv = dec.pivot_table(index="score", columns="decile", values="hit_rate")
        w("| score | " + " | ".join(str(int(c)) for c in piv.columns) + " |")
        w("|---|" + "---|" * len(piv.columns))
        for sname, row in piv.iterrows():
            w(f"| {SCORE_LABEL.get(sname, sname)} | " + " | ".join(_f(v, 2) for v in row.to_numpy()) + " |")
        w("")
    # ---------------------------------------------------------------- population adjusted
    w("## 4. Population control (big counties attract trials)")
    w("")
    w("Logistic regression of O1 (and O4) on log(1 + population within 50 km) + log(1 + population) + z(score); second "
      "model adds log(1 + prior relevant trials within 50 km). Odds ratio per SD of the score with county-bootstrap "
      "and Wald 95% CIs; likelihood-ratio test of the score; in-sample AUROC of the model without / with the score.")
    w("")
    if len(pa):
        w("| outcome | model | score | expected | n (pos) | OR per SD [bootstrap CI] | Wald CI | LR p | AUROC without -> with |")
        w("|---|---|---|---|---|---|---|---|---|")
        for _, r in pa.iterrows():
            if r.get("status") != "ok":
                w(f"| {r['outcome']} | {r['model']} | {r['score']} | | | {r.get('status')} | | | |")
                continue
            w(f"| {r['outcome']} | {r['model']} | {SCORE_LABEL.get(r['score'], r['score'])} | {r['expected_direction']} | "
              f"{int(r['n']):,} ({int(r['n_pos']):,}) | {_f(r['odds_ratio_per_sd'], 2)} [{_f(r['odds_ratio_ci_low'], 2)}, "
              f"{_f(r['odds_ratio_ci_high'], 2)}] | {_f(r['odds_ratio_wald_ci_low'], 2)}-{_f(r['odds_ratio_wald_ci_high'], 2)} | "
              f"{_f(r['lr_p'], 4)} | {_f(r['auroc_model_without_score'])} -> {_f(r['auroc_model_with_score'])} |")
    w("")
    w("Partial Spearman of each score with the post-period trial count given log population within 50 km is in "
      "`temporal_holdout_metrics.csv` (`partial_spearman_count_given_log_pop50`).")
    w("")
    # ---------------------------------------------------------------- secondary outcomes
    w("## 5. Secondary outcomes and sensitivities (T = 2021, clean variant, county): AUROC [95% CI]")
    w("")
    keys = ["composite_equal", "composite_evidence_weighted", "composite_study_partner", "composite_access_gap",
            "research_readiness_pct", "clinic_capacity_pct", "diagnostic_desert_pct", "burden_pct",
            "baseline_population_within_50km", "baseline_prior_trials_50km"]
    w("| outcome | base rate | " + " | ".join(SCORE_LABEL.get(k, k) for k in keys) + " |")
    w("|---|---|" + "---|" * len(keys))
    for o in OUTCOMES:
        s_ = sel(PRIMARY_T, "A_clean", QUERY, o)
        if not len(s_):
            continue
        w(f"| {o}: {OUTCOMES[o][1]} | {_f(s_['base_rate'].iloc[0])} | " +
          " | ".join(_ci_s(s_.loc[k], 'auroc') if k in s_.index else "NA" for k in keys) + " |")
    w("")
    w("O4 (new facilities: sites with a relevant post-T trial and no relevant pre-T trial) is the outcome that "
      "persistence of existing research sites cannot explain by itself.")
    w("")
    # ---------------------------------------------------------------- other T and variants
    w("## 6. Other cutoffs, ME/CFS alone, the not-clean variant and the leakage bound (county, O1)")
    w("")
    w("| T | variant | clean holdout | query | equal AUROC [CI] | Delta vs pop 50 km | evidence_weighted AUROC | "
      "research_readiness AUROC | desert AUROC | desert top-decile lift | population 50 km AUROC |")
    w("|---|---|---|---|---|---|---|---|---|---|---|")
    for (T, v, q), g in m[(m["outcome"] == "O1") & (m["level"] == "county")].groupby(["T", "variant", "query"],
                                                                                     sort=False):
        g = g.set_index("score")
        if "composite_equal" not in g.index:
            continue
        e = g.loc["composite_equal"]
        w(f"| {T} | {v} | {'yes' if e['clean_holdout'] else '**no (NOT a clean holdout / leakage bound)**'} | {q} | "
          f"{_ci_s(e, 'auroc')} | {_ci_s(e, 'delta_auroc_vs_population_within_50km')} | "
          f"{_ci_s(g.loc['composite_evidence_weighted'], 'auroc') if 'composite_evidence_weighted' in g.index else 'NA'} | "
          f"{_ci_s(g.loc['research_readiness_pct'], 'auroc')} | {_ci_s(g.loc['diagnostic_desert_pct'], 'auroc')} | "
          f"{_ci_s(g.loc['diagnostic_desert_pct'], 'top_decile_lift', 2)} | "
          f"{_ci_s(g.loc['baseline_population_within_50km'], 'auroc')} |")
    w("")
    a_ = sel(PRIMARY_T, "A_clean", QUERY, "O1")
    s_ = sel(PRIMARY_T, "A_static_context", QUERY, "O1")
    if len(s_):
        rows = []
        for k in ("composite_equal", "composite_evidence_weighted", "research_readiness_pct", "clinic_capacity_pct",
                  "diagnostic_desert_pct"):
            d = s_.loc[k, "auroc"] - a_.loc[k, "auroc"]
            half = (a_.loc[k, "auroc_ci_high"] - a_.loc[k, "auroc_ci_low"]) / 2
            rows.append(f"{SCORE_LABEL.get(k, k)} {d:+.3f} (CI half-width {half:.3f}; "
                        f"{'within' if abs(d) < half else 'EXCEEDS'})")
        w("**Leakage bound (plan section 2):** AUROC change when providers and facilities are NOT date-filtered: "
          + "; ".join(rows) + ".")
        w("")
    if len(lk):
        r = lk[(lk["T"] == PRIMARY_T) & (lk["variant"] == "A_clean")].iloc[0]
        w(f"What the T = {PRIMARY_T} date filters removed: {int(r['individual_npis_enumerated_after_T']):,} of "
          f"{int(r['individual_npis']):,} individual NPIs ({r['individual_npis_enumerated_after_T'] / r['individual_npis']:.1%}) "
          f"enumerated after T; {int(r['facilities_without_pre_T_date']):,} of {int(r['facilities_clinic_candidates']):,} "
          f"clinic-candidate facilities ({r['facilities_without_pre_T_date'] / r['facilities_clinic_candidates']:.1%}) "
          f"without a dated source on or before T; facility-trial rows {int(r['facility_trial_rows']):,} -> "
          f"{int(r['facility_trial_rows_pre_T']):,}; facility-NIH rows {int(r['facility_nih_rows']):,} -> "
          f"{int(r['facility_nih_rows_pre_T']):,}.")
        w("")
    # ---------------------------------------------------------------- state level
    w("## 7. State level (51 states; Spearman with post-period counts)")
    w("")
    if len(st):
        s7 = st[(st["T"] == PRIMARY_T) & (st["variant"] == "A_clean") & (st["query"] == QUERY)]
        keys7 = ["composite_equal", "composite_evidence_weighted", "research_readiness_pct", "clinic_capacity_pct",
                 "diagnostic_desert_pct", "burden_pct", "baseline_population_total"]
        w("| outcome | scale | states positive | " + " | ".join(SCORE_LABEL.get(k, k) for k in keys7) + " |")
        w("|---|---|---|" + "---|" * len(keys7))
        for (o, sc), g in s7.groupby(["outcome", "scale"], sort=False):
            g = g.set_index("score")
            w(f"| {o} | {sc} | {int(g['n_states_positive'].iloc[0])} of {int(g['n_states'].iloc[0])} | " +
              " | ".join(f"{_f(g.loc[k, 'spearman'], 2)} [{_f(g.loc[k, 'spearman_ci_low'], 2)}, "
                         f"{_f(g.loc[k, 'spearman_ci_high'], 2)}]" if k in g.index else "NA" for k in keys7) + " |")
        w("")
        w("Binary O1 at state level is not analysed: almost every state has a new relevant trial (ceiling; plan "
          "section 6).")
    w("")
    # ---------------------------------------------------------------- burden persistence
    w("## 8. Burden input persistence (CMS MMD condition 51, ME/CFS level-C proxy)")
    w("")
    if len(bp):
        w(f"| T | counties | Spearman(value at T, value {CMS_LATEST_FINAL}) [CI] | mean change (pp) | Spearman(value at T, change) | "
          "top-decile retained |")
        w("|---|---|---|---|---|---|")
        for _, r in bp.iterrows():
            w(f"| {int(r['T'])} | {int(r['n']):,} | {_f(r['spearman_cms51_T_vs_latest_final'])} [{_f(r['ci_low'])}, "
              f"{_f(r['ci_high'])}] | {_f(r['mean_change_pp'], 2)} | {_f(r['spearman_change_vs_T_value'])} | "
              f"{_f(r['top_decile_retained_share'], 2)} |")
        w("")
        w("This checks the stability of the burden input only; it is not a validation of burden (a claims proxy).")
        w("")
    # ---------------------------------------------------------------- stability
    w("## 9. Stability")
    w("")
    w("**Leave-one-state-out normalisation** (each state's counties scored against a normalisation reference of the "
      "other 50 states; plan section 7):")
    w("")
    if len(ls):
        w("| ranking | states | with >= 3 ranked counties | median within-state Spearman | min | states < 0.8 | median "
          "other-county Spearman | min | median top-25 overlap | min |")
        w("|---|---|---|---|---|---|---|---|---|---|")
        for _, r in ls.iterrows():
            below = r["states_below_0_8"] if isinstance(r["states_below_0_8"], str) and r["states_below_0_8"] else "none"
            w(f"| {r['ranking']} | {int(r['n_states'])} | {int(r['n_states_with_3plus_ranked'])} | "
              f"{_f(r['median_spearman_state_counties'])} | {_f(r['min_spearman_state_counties'])} | {below} | "
              f"{_f(r['median_spearman_other_counties'], 4)} | {_f(r['min_spearman_other_counties'], 4)} | "
              f"{_f(r['median_top25_overlap'], 1)} | {int(r['min_top25_overlap'])} |")
        w("")
    if len(lo):
        worst = lo.sort_values("spearman_state_counties").head(5)
        w("Lowest within-state agreement: " + "; ".join(
            f"{r['ranking']} {r.get('state_abbr', r['state_fips'])} ({int(r['n_state_counties_ranked'])} counties) "
            f"{_f(r['spearman_state_counties'])}" for _, r in worst.iterrows()) + ".")
        w("")
    w("**Rolling check** (equal-weight county rankings, clean variant; counties eligible in both):")
    w("")
    if len(ro):
        w("| a | b | counties | Spearman | Kendall tau-b | top-10 overlap | top-25 overlap |")
        w("|---|---|---|---|---|---|---|")
        for _, r in ro.iterrows():
            w(f"| {r['ranking_a']} | {r['ranking_b']} | {int(r['n_common']):,} | {_f(r['spearman'])} | "
              f"{_f(r['kendall_tau_b'])} | {int(r['top10_overlap'])} | {int(r['top25_overlap'])} |")
        w("")
    # ---------------------------------------------------------------- top 25
    w(f"## 10. The pre-{PRIMARY_T + 1} top-25 (equal weights) and what happened there")
    w("")
    if len(top):
        w("| pre-T rank | county | research_readiness_pct | diagnostic_desert_pct | relevant trials within 50 km <= T | "
          "new relevant trials within 50 km > T | new facilities within 50 km | new NIH core projects in county | "
          "current rank |")
        w("|---|---|---|---|---|---|---|---|---|")
        for _, r in top.iterrows():
            w(f"| {int(r['rank_equal_pre_T'])} | {r['geo_name']} | {_f(r['research_readiness_pct'], 2)} | "
              f"{_f(r['diagnostic_desert_pct'], 2)} | {int(r['prior_trials_pool_n'])} | {int(r['post_trials_pool_n'])} | "
              f"{int(r['post_new_facilities_pool_n'])} | {int(r['post_nih_new_core_in_county_n'])} | "
              f"{_f(r['current_rank_equal'], 0)} |")
        w("")
    # ---------------------------------------------------------------- limits
    w("## 11. What this cannot show, and the leakage that remains")
    w("")
    w("* **Not validated, not validatable with public data:** clinical value, diagnostic yield, patient benefit, "
      "cost, or whether a deployment in a flagged county would work. Research activity is a proxy outcome for "
      "research readiness only; trials go where investigators, funding and large populations are.")
    w("* **The two sides of the composite point in opposite directions by design.** Readiness predicts activity; the "
      "desert and vulnerability predict its absence. A composite that averages them is expected to predict activity "
      "weakly; its AUROC is not a measure of how good a deployment recommendation it is.")
    w("* **Vintage leakage kept:** SVI 2022, ACS population, RUCC 2023, county HPSA flags, NPPES taxonomies and "
      "practice addresses (current snapshot; NPIs deactivated before the snapshot are absent), the HRSA roster "
      "(active sites only), ClinicalTrials.gov site lists (a site added to a pre-T trial after T counts as pre-T), "
      "the measurement performance record (post-T data; constant within the query, it rescales expected_yield only) "
      "and the curated configs (2026). The provider / facility part is bounded by the static-context re-run "
      "(section 6); SVI / ACS cannot be bounded without earlier vintages.")
    w("* **Long COVID burden did not exist before T.** The clean variant ranks on the ME/CFS claims proxy alone; the "
      "variant with the post-period Long COVID burden is not a clean holdout.")
    w("* **NIH data start at FY2015** (pre-period NIH activity is FY2015-T); new-award status is 'first fiscal year "
      "in these data > T'.")
    w("* **Spatial dependence:** the 50 km pools overlap, so counties are not independent; the state-cluster "
      "bootstrap and within-state permutation are reported next to the county bootstrap.")
    w("")
    w("## 12. Files")
    w("")
    for k, pth in sorted(meta.get("tables", {}).items()):
        w(f"* `{pth}`")
    w(f"* `data/processed/{TABLE}.parquet` (every pre-period ranking with all components, composites and ranks, "
      "plus post-period outcome counts `post_*`)")
    if meta.get("figure"):
        w(f"* `{meta['figure']}` (draft figure)")
    w(f"* plan: `{PLAN}`; module: `src/measure_it/scoring/temporal_holdout.py`; tests: `tests/test_temporal_holdout.py`")
    w("")
    text = "\n".join(L) + "\n"
    RESULTS_MD.write_text(text)
    return str(RESULTS_MD.relative_to(config.PROJECT_ROOT))


if __name__ == "__main__":
    run()
