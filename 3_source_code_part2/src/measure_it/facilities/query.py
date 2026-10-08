"""Research-centre query for the MCP layer.

    find_relevant_research_centers(condition, measurement=None, top_n=25, include_expanded=False) -> dict

Aggregates ClinicalTrials.gov site history (`facility_trials`) and NIH RePORTER awards (`facility_nih_projects`)
per RESOLVED facility/organisation (measure_it.facilities.registry), with counts and object ids. Research
readiness is a capability signal: it is not patient burden, a trial registration is not evidence that a
measurement works, and award obligations are not expenditure. Returns UNKNOWN (with a reason) when the
condition does not resolve or has no trial or NIH record.
"""
from __future__ import annotations

from functools import lru_cache

import numpy as np
import pandas as pd

from ..config import UNKNOWN
from ..store import processed_path, read_table
from .matching import _has, resolve_condition, resolve_measurement
from .registry import SEP

FACILITY_COLUMNS = ["facility_id", "object_id", "facility_name", "primary_kind", "sources", "city", "state", "zip5",
                    "county_fips", "lat", "lon", "geocode_precision", "resolution_tier", "npi_object_ids",
                    "hrsa_site_object_ids", "member_names"]
_TABLES = ("facility_trials", "facility_nih_projects", "facilities")


def _stamp() -> tuple:
    """Modification times of the three input tables: a rebuilt table invalidates the cache."""
    return tuple(processed_path(n).stat().st_mtime_ns if processed_path(n).exists() else 0 for n in _TABLES)


@lru_cache(maxsize=2)
def _tables(stamp: tuple) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """facility_trials, facility_nih_projects and the facility columns, read once per table version (the query only
    filters them, never modifies them in place). Raises FileNotFoundError when a table is not built (not cached)."""
    return (read_table("facility_trials"), read_table("facility_nih_projects"),
            read_table("facilities", columns=FACILITY_COLUMNS))

GUARDRAILS = [
    "Research readiness (trials, NIH awards) is a capability signal, not patient burden or clinical capacity.",
    "A registered trial that mentions a measurement is not evidence that the measurement works.",
    "Award obligations are total-cost obligations summed over fiscal years (once per appl_id, no parent/subproject "
    "double counting), not expenditure, and a whole award counts when it mentions the condition.",
    "Facilities are resolved conservatively across sources; spelling variants that could not be linked stay separate, "
    "so counts per facility are lower bounds.",
    "ClinicalTrials.gov locations are ZIP or city centroids; NIH locations are the awardee organisation address.",
]


def find_relevant_research_centers(condition, measurement=None, top_n: int | None = 25,
                                   include_expanded: bool = False) -> dict:
    """Facilities with trial history or NIH awards for `condition`, aggregated per resolved facility.

    Default views: trials whose record names the condition literally (MeSH-expanded matches such as multiple
    system atrophy under dysautonomia are excluded unless include_expanded=True) and NIH projects whose
    title/abstract mention the condition (likely false positives such as pancreatic-stellate-cell "PASC" excluded;
    include_expanded=True adds terms-field-only matches). Sorted by trials + NIH core projects (a sort order for
    display, not a score); every count is returned separately.
    """
    cond = resolve_condition(condition)
    if cond["status"] != "matched":
        return {"status": UNKNOWN, "reason": f"condition: {cond['reason']}", "condition": cond, "centers": []}
    cid = cond["condition_id"]
    meas = None
    if measurement is not None and str(measurement).strip():
        meas = resolve_measurement(measurement)
        if meas["status"] != "matched":
            return {"status": UNKNOWN, "reason": f"measurement: {meas['reason']}", "condition": cond,
                    "measurement": meas, "centers": []}
    try:
        ft, fn, fac = _tables(_stamp())
    except FileNotFoundError as e:
        return {"status": UNKNOWN, "reason": f"facility tables not built: {e}", "condition": cond, "centers": []}
    tcol = "condition_ids_all" if include_expanded else "condition_ids_literal"
    ncol = "condition_ids_all" if include_expanded else "condition_ids_precision"
    t = ft[_has(ft[tcol], [cid])]
    p = fn[_has(fn[ncol], [cid])]
    if meas:
        t = t.assign(meas_hit=_has(t["measurement_classes"], meas["members"]))
    if t.empty and p.empty:
        return {"status": UNKNOWN, "reason": f"no ClinicalTrials.gov trial site or NIH RePORTER project for {cid} "
                                             f"in the resolved facility tables", "condition": cond, "centers": [],
                "guardrails": GUARDRAILS}
    from ..ingestion.nih_reporter import obligation_total

    rows = {}
    for fid, g in t.groupby("facility_id"):
        rows[fid] = {
            "n_trials": int(g["nct_id"].nunique()),
            "n_trials_open_status": int(g.loc[g["is_open_status"].astype(bool), "nct_id"].nunique()),
            "n_trials_site_recruiting": int(g.loc[g["site_recruiting"].astype(bool), "nct_id"].nunique()),
            "n_interventional": int(g.loc[g["study_type"] == "INTERVENTIONAL", "nct_id"].nunique()),
            "n_device_or_diagnostic_trials": int(g.loc[g["is_device_or_diagnostic"].astype(bool), "nct_id"].nunique()),
            "latest_start_year": (int(g["start_year"].max()) if g["start_year"].notna().any() else None),
            "trial_object_ids": sorted("trial:" + g["nct_id"].unique()),
            "measurement_classes_registered": sorted({c for s in g["measurement_classes"] for c in str(s).split(SEP)
                                                      if c}),
        }
        if meas:
            rows[fid]["n_trials_mentioning_measurement"] = int(g.loc[g["meas_hit"], "nct_id"].nunique())
    for fid, g in p.groupby("facility_id"):
        pis = {i.strip() for s in g["pi_profile_ids"].dropna() for i in str(s).split(";") if i.strip()}
        r = rows.setdefault(fid, {})
        r.update({
            "n_nih_projects": int(g["appl_id"].nunique()),
            "n_nih_core_projects": int(g["core_project_num"].nunique()),
            "n_nih_active_core_projects": int(g.loc[g["is_active"].fillna(False).astype(bool),
                                                    "core_project_num"].nunique()),
            "nih_award_obligations_usd": (float(obligation_total(g)) if g["award_amount"].notna().any() else None),
            "n_nih_award_amount_missing": int(g["award_amount"].isna().sum()),
            "n_nih_pis": len(pis),
            "nih_fiscal_years": f"{int(g['fiscal_year'].min())}-{int(g['fiscal_year'].max())}",
            "nih_org_names": sorted(g["org_name"].dropna().unique().tolist()),
            "nih_object_ids": sorted("nih:" + g["appl_id"].astype(str).unique()),
        })
    df = pd.DataFrame.from_dict(rows, orient="index")
    df.index.name = "facility_id"
    df = df.reset_index().merge(fac, on="facility_id", how="left")
    for c in ("n_trials", "n_trials_open_status", "n_trials_site_recruiting", "n_interventional",
              "n_device_or_diagnostic_trials", "n_nih_projects", "n_nih_core_projects", "n_nih_active_core_projects",
              "n_nih_pis", "n_nih_award_amount_missing"):
        if c not in df:
            df[c] = 0
        df[c] = df[c].fillna(0).astype(int)
    if meas:
        df["n_trials_mentioning_measurement"] = df.get("n_trials_mentioning_measurement", 0)
        df["n_trials_mentioning_measurement"] = df["n_trials_mentioning_measurement"].fillna(0).astype(int)
    for c in ("trial_object_ids", "nih_object_ids", "measurement_classes_registered", "nih_org_names"):
        if c not in df:
            df[c] = [[] for _ in range(len(df))]
        df[c] = df[c].map(lambda v: v if isinstance(v, list) else [])
    df["_sort"] = df["n_trials"] + df["n_nih_core_projects"]
    sort_cols = ["_sort", "n_trials", "n_nih_core_projects", "facility_id"]
    if meas:
        sort_cols = ["n_trials_mentioning_measurement"] + sort_cols
    df = df.sort_values(sort_cols, ascending=[False] * (len(sort_cols) - 1) + [True]).drop(columns="_sort")
    n_total = len(df)
    if top_n:
        df = df.head(int(top_n))
    centers = []
    for r in df.itertuples(index=False):
        d = r._asdict()
        for k, v in list(d.items()):
            if isinstance(v, (np.integer,)):
                d[k] = int(v)
            elif isinstance(v, (np.floating,)):
                d[k] = None if np.isnan(v) else float(v)
            elif v is pd.NA:
                d[k] = None
        d["npi_object_ids"] = [x for x in str(d.get("npi_object_ids") or "").split(SEP) if x][:10]
        d["hrsa_site_object_ids"] = [x for x in str(d.get("hrsa_site_object_ids") or "").split(SEP) if x]
        centers.append(d)
    return {
        "status": "ok", "condition": cond, "measurement": meas,
        "view": ("MeSH/terms-expanded" if include_expanded else
                 "literal trial condition match + NIH title/abstract match (precision view)"),
        "n_facilities": n_total, "n_returned": len(centers),
        "n_trials_total": int(t["nct_id"].nunique()), "n_nih_projects_total": int(p["appl_id"].nunique()),
        "centers": centers, "guardrails": GUARDRAILS,
        "provenance": {"tables": ["facility_trials", "facility_nih_projects", "facilities"],
                       "producer": "measure_it.facilities.registry"},
    }


def research_centers_frame(result: dict) -> pd.DataFrame:
    cols = ["facility_id", "object_id", "facility_name", "primary_kind", "city", "state", "n_trials",
            "n_trials_open_status", "n_device_or_diagnostic_trials", "latest_start_year", "n_nih_projects",
            "n_nih_core_projects", "n_nih_active_core_projects", "nih_award_obligations_usd", "n_nih_pis",
            "geocode_precision", "resolution_tier"]
    df = pd.DataFrame(result.get("centers", []))
    return df[[c for c in cols if c in df.columns]] if len(df) else pd.DataFrame(columns=cols)


# --------------------------------------------------------------------------------------------------------------
# example outputs for the write-up (results/tables/facility_*_examples.csv, draft overlay map)
# --------------------------------------------------------------------------------------------------------------

EXAMPLE_QUERIES = [  # fixed before looking at any output: SPEC example county, a large metro, a rural county
    ("06073", "Long COVID", "wearable autonomic monitoring"),
    ("17031", "ME/CFS", "wearable autonomic monitoring"),
    ("21095", "Long COVID", "wearable autonomic monitoring"),
    ("06073", "POTS", "autonomic testing"),
]


def write_examples() -> None:
    from ..config import FIGURES, TABLES
    from .matching import candidates_frame, find_candidate_clinics

    TABLES.mkdir(parents=True, exist_ok=True)
    frames = []
    for geo, cond, meas in EXAMPLE_QUERIES:
        r = find_candidate_clinics(geo, cond, meas)
        df = candidates_frame(r)
        if df.empty:
            df = pd.DataFrame([{"status": r["status"], "reason": r.get("reason")}])
        df.insert(0, "query_measurement", meas)
        df.insert(0, "query_condition", cond)
        df.insert(0, "query_geography", geo)
        df.insert(3, "status", r["status"])
        df.insert(4, "n_facilities_in_pool", r.get("n_facilities_in_pool"))
        df.insert(5, "n_eligible", r.get("n_eligible"))
        for c, v in (r.get("n_with_characteristic") or {}).items():
            df[f"pool_n_with__{c}"] = v
        df["absence_note"] = r.get("absence_note", "")
        frames.append(df)
    pd.concat(frames, ignore_index=True).to_csv(TABLES / "facility_candidate_examples.csv", index=False)
    rows = []
    for cond in ("long_covid", "me_cfs", "pots", "dysautonomia"):
        for expanded in (False, True):
            r = find_relevant_research_centers(cond, top_n=None, include_expanded=expanded)
            df = research_centers_frame(r)
            df.insert(0, "view", "expanded" if expanded else "precision")
            df.insert(0, "condition_id", cond)
            df["rank_in_display_order"] = np.arange(1, len(df) + 1)
            rows.append(df)
    allc = pd.concat(rows, ignore_index=True)
    allc[allc["rank_in_display_order"] <= 25].to_csv(TABLES / "facility_research_centers_examples.csv", index=False)
    summ = allc.groupby(["condition_id", "view"]).agg(
        n_facilities=("facility_id", "size"), n_with_trials=("n_trials", lambda x: int((x > 0).sum())),
        n_with_nih=("n_nih_projects", lambda x: int((x > 0).sum())),
        n_with_both=("facility_id", lambda x: int(((allc.loc[x.index, "n_trials"] > 0) &
                                                   (allc.loc[x.index, "n_nih_projects"] > 0)).sum())),
        n_states=("state", "nunique")).reset_index()
    summ.to_csv(TABLES / "facility_research_centers_summary.csv", index=False)
    _overlay_map(FIGURES / "drafts" / "facility_infrastructure_overlay.png")


def _overlay_map(path) -> None:
    import geopandas as gpd
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from ..config import PROCESSED
    st = gpd.read_parquet(PROCESSED / "state_boundaries_2024.geoparquet")
    st = st[~st["STATEFP"].isin(["02", "15", "60", "66", "69", "72", "78"])] if "STATEFP" in st.columns else st
    f = read_table("facilities", columns=["facility_id", "lat", "lon", "is_hrsa_site", "sources"])
    rr = read_table("research_site_registry", columns=["facility_id", "lat", "lon", "n_trials_literal",
                                                       "n_nih_projects_precision"])
    box = dict(lon=(-125, -66), lat=(24, 50))

    def inbox(d):
        return d[d["lon"].between(*box["lon"]) & d["lat"].between(*box["lat"])]
    fig, ax = plt.subplots(figsize=(13, 7.5))
    st.boundary.plot(ax=ax, color="#bbbbbb", linewidth=0.4)
    h = inbox(f[f["is_hrsa_site"].astype(bool)])
    ax.scatter(h["lon"], h["lat"], s=1.5, c="#3a9e5c", alpha=0.35, label=f"HRSA health-center sites (n={len(h):,})")
    t = inbox(rr[rr["n_trials_literal"] > 0])
    ax.scatter(t["lon"], t["lat"], s=3 + 2 * np.sqrt(t["n_trials_literal"]), facecolors="none", edgecolors="#1f5fa8",
               linewidths=0.5, alpha=0.6, label=f"trial sites, target conditions (n={len(t):,}; size ~ trials)")
    n = inbox(rr[rr["n_nih_projects_precision"] > 0])
    ax.scatter(n["lon"], n["lat"], s=12, marker="^", c="#c0392b", alpha=0.7,
               label=f"NIH RePORTER awardees, target conditions (n={len(n):,})")
    ax.set_xlim(*box["lon"])
    ax.set_ylim(*box["lat"])
    ax.set_axis_off()
    ax.legend(loc="lower left", fontsize=8, frameon=False)
    ax.set_title("Draft: facility and research infrastructure (contiguous U.S.; resolved facilities; "
                 "ClinicalTrials.gov locations are ZIP/city centroids)", fontsize=10)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    write_examples()
