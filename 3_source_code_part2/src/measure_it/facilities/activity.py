"""Measurement-specific clinician capacity and experience sites (stage F, docs/HARMONIZATION_CONTRACT.md section 5).

Why. The scoring's ``clinic_capacity`` counts NPPES providers/facilities whose SELF-REPORTED taxonomy is in a
measurement's implementer groups (configs/relevance.yaml). A location shuffle leaves it almost unchanged (Spearman
0.965, SCORING_RESULTS Test 6b): it measures provider density, not who performs the measurement. This module
replaces "could plausibly" with "billed it": Medicare FFS claims for the measurement's HCPCS/CPT codes
(``provider_measurement_activity``, configs/measurement_hcpcs.yaml) joined to NPPES ``providers`` for county and
specialty.

Generic over measurement ids: every class of configs/measurements.yaml and every bundle of configs/relevance.yaml
(``measurement_bundles``; a bundle is the union of its member classes). A new measurement needs config rows only.

Outputs
-------
``geo_measurement_capacity``  geo (county / state) x measurement id: active clinicians (individual rendering NPIs
    with >= 1 claim row of the measurement's capacity basis), organisations, services, beneficiary sums, per 100k
    adults, county-level within-50-km counts and ``measurement_capacity`` / ``_pct`` (percentile composite of the
    clinician rates, the same pr() as the scoring). ``capacity_basis`` = ``dedicated`` when any member class has a
    dedicated code in the config, else ``related``, else ``proxy`` (labelled; e.g. nailfold capillaroscopy has no
    CPT code), else ``none`` (all-NaN capacity: no claims signal exists for that measurement).
``measurement_experience_sites``  U.S. facilities with a registered trial (facility_trials.measurement_classes) or an
    NIH project (measurement_mentions on the project's text) that mentions a member class; with NCT / appl ids and
    NIH contact PIs. A registered study is research experience, not evidence the site offers the test clinically.
``provider_measurement_summary``  NPI x measurement id: best code role, services per role, codes billed, county.

Reports (results/tables): ``stage_f_capacity_vs_clinic_capacity.csv`` (Spearman of the new capacity with the scoring's
``clinic_capacity`` per bundle) and ``stage_f_capacity_location_shuffle.csv`` (Test-6b analogue: individual-NPI
locations permuted nationally / within state; Spearman of the shuffled capacity with the observed one).

Run: ``uv run python -m measure_it.facilities.activity [--draws 200 --state-draws 100]`` (``build()`` = tables only,
``report()`` = comparison + shuffle, which reads deployment_opportunities).
"""
from __future__ import annotations

import argparse
import json
import time
from functools import lru_cache

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import spearmanr

from ..config import SEED, TABLES, load_config, utc_now_iso
from ..geography.crosswalk import STATE_FIPS, zip_to_geo
from ..provenance import add_provenance
from ..store import read_table, table_exists, write_table
from .registry import haversine_km

PRODUCER = "facilities.activity"
SEP = "|"
ROLES = ("dedicated", "related", "proxy")
LEVELS = ("county", "state")
SUMMARY_TABLE = "provider_measurement_summary"
CAPACITY_TABLE = "geo_measurement_capacity"
SITES_TABLE = "measurement_experience_sites"

CAPACITY_NOTE = (
    "Clinicians = individual rendering NPIs (CMS entity I) billing >= 1 code of the measurement's capacity basis to "
    "Medicare FFS Part B in the data year, located by their NPPES primary practice ZIP (ZCTA internal point -> 2024 "
    "county; CMS ZIP when the NPI is not in providers). Medicare FFS only; CMS suppresses provider x code rows with "
    "< 11 beneficiaries, so counts are lower bounds. Billing a code does not mean seeing Long COVID/ME/CFS or any "
    "target condition. capacity_basis 'proxy' = no dedicated code exists; the count is an indirect signal. Per-100k "
    "rates in the geography use ACS adults 18+; within-50-km rates use the total population within 50 km."
)


# ------------------------------------------------------------------------------------------------------------------
# measurement ids
# ------------------------------------------------------------------------------------------------------------------

@lru_cache(maxsize=1)
def measurement_members() -> dict[str, dict]:
    """measurement id -> {kind, members}: every class of measurements.yaml and every bundle of relevance.yaml."""
    out = {c["id"]: {"kind": "class", "members": [c["id"]]}
           for c in load_config("measurements")["measurement_classes"]}
    for bid, b in (load_config("relevance").get("measurement_bundles") or {}).items():
        out[bid] = {"kind": "bundle", "members": list(b["members"])}
    return out


def members_of(measurement_id: str) -> list[str]:
    mm = measurement_members()
    if measurement_id not in mm:
        raise KeyError(f"unknown measurement id {measurement_id!r} (configs/measurements.yaml, relevance.yaml)")
    return mm[measurement_id]["members"]


@lru_cache(maxsize=1)
def configured_roles() -> dict[str, set]:
    """measurement class -> set of code roles configured for it (dedicated / related / proxy)."""
    from ..ingestion.cms_physician_service import code_map, load_codes_config
    cm = code_map(load_codes_config())
    return {k: set(v) for k, v in cm.groupby("measurement_class")["code_role"]}


def capacity_basis(measurement_id: str, roles: dict | None = None) -> str:
    """dedicated > related > proxy > none, over the member classes' CONFIGURED codes (not over what was billed)."""
    roles = configured_roles() if roles is None else roles
    have = set().union(*[roles.get(m, set()) for m in members_of(measurement_id)])
    for r in ROLES:
        if r in have:
            return r
    return "none"


# ------------------------------------------------------------------------------------------------------------------
# provider-level summary
# ------------------------------------------------------------------------------------------------------------------

@lru_cache(maxsize=1)
def activity() -> pd.DataFrame:
    return read_table("provider_measurement_activity")


@lru_cache(maxsize=1)
def provider_locations() -> pd.DataFrame:
    """Every NPI in provider_measurement_activity with a location: NPPES providers first, CMS ZIP5 as fallback."""
    a = activity()
    first = a.drop_duplicates("npi").set_index("npi")[["entity_code", "cms_zip5", "cms_state_fips", "cms_provider_type"]]
    p = read_table("providers", columns=["npi", "entity_type", "county_fips", "state_fips", "lat", "lon",
                                         "specialty_groups", "primary_taxonomy_grouping", "deactivation_date"])
    p = p[p["npi"].isin(first.index)].set_index("npi")
    loc = first.join(p, how="left")
    loc["location_source"] = np.where(loc["lat"].notna(), "nppes_providers", None)
    need = loc["lat"].isna()
    if need.any():
        z = zip_to_geo(loc.loc[need, "cms_zip5"].fillna(""))
        loc.loc[need, "lat"] = z["lat"].values
        loc.loc[need, "lon"] = z["lon"].values
        loc.loc[need, "county_fips"] = z["county_fips"].values
        loc.loc[need, "state_fips"] = loc.loc[need, "cms_state_fips"].values
        loc.loc[need & loc["lat"].notna(), "location_source"] = "cms_zip5"
    loc["in_providers"] = loc["entity_type"].notna()
    return loc.reset_index()


def provider_summary(measurement_id: str) -> pd.DataFrame:
    """NPI x measurement id: services per role, best role, codes billed (code:services), location."""
    a = activity()
    a = a[a["measurement_class"].isin(members_of(measurement_id))]
    if a.empty:
        return pd.DataFrame(columns=["npi", "measurement_id"])
    a = a.drop_duplicates(["npi", "hcpcs_cd", "code_role"])  # a code shared by two member classes counts once
    a = a.sort_values("code_role", key=lambda s: s.map({r: i for i, r in enumerate(ROLES)}))
    a = a.drop_duplicates(["npi", "hcpcs_cd"])
    piv = a.pivot_table(index="npi", columns="code_role", values="tot_services", aggfunc="sum", fill_value=0)
    out = pd.DataFrame(index=piv.index)
    for r in ROLES:
        out[f"services_{r}"] = piv[r] if r in piv else 0.0
    out["best_code_role"] = np.select([out[f"services_{r}"] > 0 for r in ROLES], list(ROLES), default="none")
    g = a.groupby("npi")
    out["tot_services"] = g["tot_services"].sum()
    out["tot_benes_sum"] = g["tot_benes_sum_over_pos"].sum()
    out["n_codes"] = g["hcpcs_cd"].nunique()
    a = a.assign(_lab=a["hcpcs_cd"] + ":" + a["tot_services"].round().astype(int).astype(str)
                 + np.where(a["code_role"] == "proxy", " (proxy)", np.where(a["code_role"] == "related",
                                                                           " (related)", "")))
    out["codes_billed"] = a.sort_values("tot_services", ascending=False).groupby("npi")["_lab"].agg("; ".join)
    out["entity_code"] = g["entity_code"].first()
    out["cms_provider_type"] = g["cms_provider_type"].first()
    out["data_year"] = g["data_year"].first()
    out = out.reset_index()
    loc = provider_locations()[["npi", "county_fips", "state_fips", "lat", "lon", "location_source", "in_providers"]]
    out = out.merge(loc, on="npi", how="left")
    out.insert(1, "measurement_id", measurement_id)
    return out


# ------------------------------------------------------------------------------------------------------------------
# geography
# ------------------------------------------------------------------------------------------------------------------

def percentile_rank(x) -> np.ndarray:
    """pr(x) = rank(average ties) / n over non-null values (the scoring's and geography module's definition)."""
    return pd.Series(np.asarray(x, dtype=float)).rank(pct=True, method="average").to_numpy()


@lru_cache(maxsize=2)
def universe(level: str) -> pd.DataFrame:
    cols = ["geo_id", "geo_level", "geo_name", "state_fips", "state_abbr", "lat", "lon", "population_total",
            "population_adults_18plus", "population_within_50km", "small_population_flag", "in_analysis_universe",
            "condition_id"]
    f = read_table("geo_condition_features", columns=cols)
    f = f[(f["condition_id"] == "long_covid") & (f["geo_level"] == level) & f["in_analysis_universe"]]
    return f.drop(columns=["condition_id"]).sort_values("geo_id").reset_index(drop=True)


def radius_km() -> float:
    return float(load_config("scoring").get("clinic_matching", {}).get("default_radius_km", 50))


@lru_cache(maxsize=1)
def location_slots() -> tuple[pd.DataFrame, np.ndarray, pd.Series]:
    """Location slots for the shuffle: every geocoded individual NPI of `providers` plus activity NPIs located only by
    their CMS ZIP. Returns (slot table [lat, lon, county_fips, state_fips], npi -> slot index array, npi index)."""
    p = read_table("providers", columns=["npi", "entity_type", "county_fips", "state_fips", "lat", "lon"])
    p = p[(p["entity_type"].astype(str) == "1") & p["lat"].notna()][["npi", "lat", "lon", "county_fips", "state_fips"]]
    extra = provider_locations()
    extra = extra[(extra["location_source"] == "cms_zip5") & (extra["entity_code"] == "I")]
    p = pd.concat([p, extra[["npi", "lat", "lon", "county_fips", "state_fips"]]], ignore_index=True)
    p = p.drop_duplicates("npi").reset_index(drop=True)
    for c in ("county_fips", "state_fips"):
        p[c] = p[c].astype(object).where(p[c].notna(), "")
    key = pd.MultiIndex.from_frame(p[["lat", "lon", "county_fips", "state_fips"]].astype(
        {"lat": float, "lon": float}).round({"lat": 6, "lon": 6}))
    codes, uniq = pd.factorize(key)
    slots = uniq.to_frame(index=False)
    slots.columns = ["lat", "lon", "county_fips", "state_fips"]
    return slots, codes.astype(np.int64), pd.Series(np.arange(len(p)), index=p["npi"].astype(str))


def _incidence(rows, cols, shape) -> sparse.csr_matrix:
    m = sparse.csr_matrix((np.ones(len(rows), dtype=np.float32), (rows, cols)), shape=shape)
    m.sum_duplicates()
    m.data[:] = 1.0
    return m


@lru_cache(maxsize=2)
def geo_incidence(level: str) -> tuple[sparse.csr_matrix, sparse.csr_matrix | None]:
    """(slot-in-geo, slot-within-radius) incidence matrices geo x slot."""
    slots, _, _ = location_slots()
    geo = universe(level)
    key = "county_fips" if level == "county" else "state_fips"
    gi = pd.Series(np.arange(len(geo)), index=geo["geo_id"])
    r = slots[key].map(gi)
    ok = r.notna().to_numpy()
    loc_in = _incidence(r[ok].to_numpy(int), np.flatnonzero(ok), (len(geo), len(slots)))
    if level != "county":
        return loc_in, None
    from sklearn.neighbors import BallTree
    rad = radius_km()
    tree = BallTree(np.radians(slots[["lat", "lon"]].to_numpy(float)), metric="haversine")
    hits = tree.query_radius(np.radians(geo[["lat", "lon"]].to_numpy(float)), r=(rad + 0.5) / 6371.0088)
    rows, cols = [], []
    for i, h in enumerate(hits):
        if len(h):
            d = haversine_km(geo["lat"].iloc[i], geo["lon"].iloc[i], slots["lat"].to_numpy()[h],
                             slots["lon"].to_numpy()[h])
            h = h[d <= rad]
            rows.append(np.full(len(h), i))
            cols.append(h)
    loc_50 = _incidence(np.concatenate(rows), np.concatenate(cols), (len(geo), len(slots)))
    return loc_in, loc_50


def _rate(n, den) -> np.ndarray:
    den = np.asarray(den, dtype=float)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(den > 0, np.asarray(n, float) / den * 1e5, np.nan)


def capacity_counts(level: str, npi_weight: dict[str, np.ndarray], npi_slot: np.ndarray) -> dict[str, np.ndarray]:
    """Aggregate per-NPI weights (aligned with the slot NPI index) to geographies: in-geo and (county) within radius."""
    loc_in, loc_50 = geo_incidence(level)
    n_slots = loc_in.shape[1]
    out = {}
    for name, w in npi_weight.items():
        per_slot = np.bincount(npi_slot, weights=w, minlength=n_slots)
        out[f"{name}_in_geo"] = np.asarray(loc_in @ per_slot).ravel()
        if loc_50 is not None:
            out[f"{name}_within_50km"] = np.asarray(loc_50 @ per_slot).ravel()
    return out


def capacity_components(level: str, cnt: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    geo = universe(level)
    out = {"active_clinicians_per_100k_adults": _rate(cnt["clin_in_geo"], geo["population_adults_18plus"])}
    rates = [out["active_clinicians_per_100k_adults"]]
    if "clin_within_50km" in cnt:
        out["active_clinicians_within_50km_per_100k"] = _rate(cnt["clin_within_50km"], geo["population_within_50km"])
        rates.append(out["active_clinicians_within_50km_per_100k"])
    pcts = [percentile_rank(r) for r in rates]
    with np.errstate(invalid="ignore"):
        cap = np.nanmean(np.vstack(pcts), axis=0) if not all(np.isnan(p).all() for p in pcts) else \
            np.full(len(geo), np.nan)
    out["measurement_capacity"] = cap
    out["measurement_capacity_pct"] = percentile_rank(cap)
    return out


def _npi_weights(summary: pd.DataFrame, basis: str) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    """Per-slot-NPI weight vectors for the basis tier (individual clinicians) plus organisations/services."""
    _, slot_codes, npi_index = location_slots()
    n = len(npi_index)
    clin = {r: np.zeros(n) for r in ROLES}
    serv = np.zeros(n)
    benes = np.zeros(n)
    ind = summary[(summary["entity_code"] == "I") & summary["npi"].isin(npi_index.index)]
    idx = npi_index.loc[ind["npi"]].to_numpy()
    for r in ROLES:
        clin[r][idx] = (ind[f"services_{r}"].to_numpy() > 0).astype(float)
    if basis in ROLES:
        serv[idx] = ind[f"services_{basis}"].to_numpy(float)
        benes[idx] = np.where(ind[f"services_{basis}"].to_numpy() > 0, ind["tot_benes_sum"].fillna(0).to_numpy(float), 0)
    w = {"clin": clin[basis] if basis in ROLES else np.zeros(n), "serv": serv, "benes": benes}
    w.update({f"clin_{r}": clin[r] for r in ROLES})
    return w, {"slot": slot_codes}


def _org_counts(level: str, summary: pd.DataFrame, basis: str) -> np.ndarray:
    geo = universe(level)
    key = "county_fips" if level == "county" else "state_fips"
    if basis not in ROLES:
        return np.zeros(len(geo))
    o = summary[(summary["entity_code"] == "O") & (summary[f"services_{basis}"] > 0)]
    c = o.groupby(key)["npi"].nunique()
    return c.reindex(geo["geo_id"]).fillna(0).to_numpy(float)


def build_capacity(measurement_ids: list[str] | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    ids = measurement_ids or list(measurement_members())
    _, slot_codes, _ = location_slots()
    frames, summaries = [], []
    for mid in ids:
        basis = capacity_basis(mid)
        summ = provider_summary(mid)
        if len(summ):
            summaries.append(summ)
        else:
            summ = pd.DataFrame(columns=["npi", "entity_code", "county_fips", "state_fips", "tot_benes_sum"]
                                + [f"services_{r}" for r in ROLES])
        w, _ = _npi_weights(summ, basis)
        for level in LEVELS:
            geo = universe(level)
            cnt = capacity_counts(level, w, slot_codes)
            comp = capacity_components(level, cnt)
            if basis == "none":
                comp = {k: np.full(len(geo), np.nan) for k in comp}
            df = geo[["geo_id", "geo_level", "geo_name", "state_abbr", "population_total",
                      "population_adults_18plus", "small_population_flag"]].copy()
            df["measurement_id"] = mid
            df["measurement_kind"] = measurement_members()[mid]["kind"]
            df["member_classes"] = SEP.join(members_of(mid))
            df["capacity_basis"] = basis
            df["active_clinicians_n"] = cnt["clin_in_geo"]
            for r in ROLES:
                df[f"clinicians_{r}_n"] = cnt[f"clin_{r}_in_geo"]
            df["active_organizations_n"] = _org_counts(level, summ, basis)
            df["services_n"] = cnt["serv_in_geo"]
            df["beneficiaries_sum"] = cnt["benes_in_geo"]
            df["services_per_100k_adults"] = _rate(cnt["serv_in_geo"], geo["population_adults_18plus"])
            if level == "county":
                df["active_clinicians_within_50km_n"] = cnt["clin_within_50km"]
                df["population_within_50km"] = geo["population_within_50km"].to_numpy()
            else:
                df["active_clinicians_within_50km_n"] = np.nan
                df["population_within_50km"] = np.nan
            for k, v in comp.items():
                df[k] = v
            if "active_clinicians_within_50km_per_100k" not in df:
                df["active_clinicians_within_50km_per_100k"] = np.nan
            frames.append(df)
    cap = pd.concat(frames, ignore_index=True)
    year = int(activity()["data_year"].iloc[0])
    cap["data_year"] = year
    cap["object_id"] = "geo:" + cap["geo_id"]
    for c in ["active_clinicians_n", "active_organizations_n", "active_clinicians_within_50km_n"] + \
            [f"clinicians_{r}_n" for r in ROLES]:
        cap[c] = cap[c].round().astype("Int64")
    cap = add_provenance(
        cap, data_layer="derived",
        source_name="CMS Medicare Physician & Other Practitioners by Provider and Service + NPPES providers + ACS",
        source_version=f"Medicare data year {year}; configs/measurement_hcpcs.yaml",
        retrieved_at=utc_now_iso(), evidence_type="derived_score",
        source_record_id=lambda d: d["measurement_id"] + "|" + d["geo_id"],
        source_geographic_resolution="county", evidence_level="administrative_claims_ffs",
        provenance_notes=CAPACITY_NOTE)
    cap.loc[cap["geo_level"] == "state", "source_geographic_resolution"] = "state"
    summ = pd.concat(summaries, ignore_index=True) if summaries else pd.DataFrame()
    if len(summ):
        summ["object_id"] = "npi:" + summ["npi"]
        summ = add_provenance(
            summ, data_layer="facility", source_name="CMS Medicare Physician & Other Practitioners by Provider and "
            "Service (provider_measurement_activity) + NPPES providers", source_version=f"Medicare data year {year}",
            retrieved_at=utc_now_iso(), evidence_type="administrative_claims",
            source_record_id=lambda d: d["npi"] + "|" + d["measurement_id"],
            source_geographic_resolution="zcta_centroid", evidence_level="administrative_claims_ffs",
            provenance_notes=CAPACITY_NOTE)
    return cap, summ


# ------------------------------------------------------------------------------------------------------------------
# experience sites
# ------------------------------------------------------------------------------------------------------------------

def _has(s: pd.Series, ids) -> pd.Series:
    ids = set(ids)
    return s.fillna("").map(lambda v: bool(ids & set(str(v).split(SEP))))


@lru_cache(maxsize=1)
def _nih_mentions() -> pd.DataFrame:
    mm = read_table("measurement_mentions", columns=["source", "record_id", "measurement_id", "evidence_flag",
                                                     "objective_flag", "nonhuman_context_hit"])
    mm = mm[(mm["source"] == "nih") & mm["evidence_flag"] & mm["objective_flag"]]
    return mm.groupby("record_id")["measurement_id"].agg(lambda s: SEP.join(sorted(set(s)))).rename(
        "nih_measurement_classes").reset_index()


def experience_sites(measurement_id: str) -> pd.DataFrame:
    """U.S. facilities with >= 1 registered trial or NIH project whose text mentions a member class."""
    members = members_of(measurement_id)
    ft = read_table("facility_trials", columns=["facility_id", "nct_id", "measurement_classes",
                                                "condition_ids_literal", "start_year", "overall_status",
                                                "brief_title"])
    ft = ft[_has(ft["measurement_classes"], members)]
    fn = read_table("facility_nih_projects", columns=["facility_id", "appl_id", "core_project_num", "fiscal_year",
                                                      "contact_pi_name", "project_title", "condition_ids_precision"])
    fn["appl_id"] = fn["appl_id"].astype(str)
    fn = fn.merge(_nih_mentions(), left_on="appl_id", right_on="record_id", how="inner")
    fn = fn[_has(fn["nih_measurement_classes"], members)]
    fac = read_table("facilities", columns=["facility_id", "object_id", "facility_name", "primary_kind", "address",
                                            "city", "state", "zip5", "county_fips", "state_fips", "lat", "lon",
                                            "geocode_precision"])
    ids = sorted(set(ft["facility_id"]) | set(fn["facility_id"]))
    # U.S. only: a facility with a U.S. state FIPS (the registries hold U.S. records; territories kept)
    f = fac[fac["facility_id"].isin(ids) & fac["state_fips"].isin(list(STATE_FIPS))].copy()

    def matched(classes: pd.Series) -> str:
        return SEP.join(sorted({c for v in classes for c in str(v).split(SEP)} & set(members)))

    t = ft.groupby("facility_id").agg(
        n_trials=("nct_id", "nunique"), nct_ids=("nct_id", lambda s: SEP.join(sorted(set(s)))),
        trial_classes=("measurement_classes", matched),
        trial_conditions=("condition_ids_literal", lambda s: SEP.join(sorted({c for v in s for c in str(v).split(SEP) if c}))),
        latest_trial_start_year=("start_year", "max"))
    n = fn.groupby("facility_id").agg(
        n_nih_projects=("appl_id", "nunique"), n_nih_core_projects=("core_project_num", "nunique"),
        appl_ids=("appl_id", lambda s: SEP.join(sorted(set(s)))),
        nih_classes=("nih_measurement_classes", matched),
        nih_contact_pis=("contact_pi_name", lambda s: SEP.join(sorted({str(x) for x in s if pd.notna(x)}))),
        nih_conditions=("condition_ids_precision", lambda s: SEP.join(sorted({c for v in s for c in str(v).split(SEP) if c}))),
        latest_nih_fiscal_year=("fiscal_year", "max"))
    f = f.merge(t, on="facility_id", how="left").merge(n, on="facility_id", how="left")
    for c in ("n_trials", "n_nih_projects", "n_nih_core_projects"):
        f[c] = f[c].fillna(0).astype(int)
    f["member_classes_matched"] = [SEP.join(sorted({x for x in (str(a) + SEP + str(b)).split(SEP)
                                                    if x and x != "nan"})) for a, b in zip(f["trial_classes"],
                                                                                         f["nih_classes"])]
    f["condition_ids"] = [SEP.join(sorted({x for x in (str(a) + SEP + str(b)).split(SEP) if x and x != "nan"}))
                          for a, b in zip(f["trial_conditions"], f["nih_conditions"])]
    f.insert(0, "measurement_id", measurement_id)
    f = f.rename(columns={"object_id": "facility_object_id"})
    f["evidence_basis"] = np.select([(f["n_trials"] > 0) & (f["n_nih_projects"] > 0), f["n_trials"] > 0],
                                    ["trial+nih", "trial"], default="nih")
    return f.drop(columns=["trial_classes", "nih_classes", "trial_conditions", "nih_conditions"])


def build_experience_sites(measurement_ids: list[str] | None = None) -> pd.DataFrame:
    ids = measurement_ids or list(measurement_members())
    frames = [experience_sites(m) for m in ids]
    df = pd.concat([f for f in frames if len(f)], ignore_index=True)
    df["object_id"] = df["facility_object_id"]
    return add_provenance(
        df, data_layer="facility", source_name="facility_trials (ClinicalTrials.gov) + facility_nih_projects "
        "(NIH RePORTER) + measurement_mentions", source_version="processed tables of this repository",
        retrieved_at=utc_now_iso(), evidence_type="text_mined",
        source_record_id=lambda d: d["measurement_id"] + "|" + d["facility_id"],
        source_geographic_resolution="zcta_centroid", evidence_level="registered_study_mention",
        provenance_notes=("A registered trial or NIH project whose text mentions the measurement class at this "
                          "facility (U.S. only): research experience, not evidence that the facility offers the test "
                          "clinically, nor that the measurement works. Facility-record linkage limits apply "
                          "(docs/FACILITY_MATCHING.md)."))


# ------------------------------------------------------------------------------------------------------------------
# reports: comparison with clinic_capacity; location shuffle
# ------------------------------------------------------------------------------------------------------------------

def compare_with_clinic_capacity(cap: pd.DataFrame) -> pd.DataFrame:
    if not table_exists("deployment_opportunities"):
        return pd.DataFrame()
    from ..config import PROCESSED
    cols = ["condition_id", "measurement_id", "geo_level", "geo_id", "clinic_capacity", "clinic_capacity_pct",
            "cc_impl_providers_in_geo_per_100k", "rank_eligible"]
    have = set(json.loads((PROCESSED / "deployment_opportunities.meta.json").read_text()).get("columns", {}))
    density = "clinic_capacity_density" in have
    if density:  # since 2026-10-07 clinic_capacity may BE the activity measure: compare with the density version
        cols += ["clinic_capacity_density"]
    d = read_table("deployment_opportunities", columns=cols)
    if density:
        d["clinic_capacity"] = d["clinic_capacity_density"]
    d = d[d["condition_id"] == "long_covid"]  # clinic_capacity does not depend on the condition
    rows = []
    for (mid, level), g in d.groupby(["measurement_id", "geo_level"]):
        c = cap[(cap["measurement_id"] == mid) & (cap["geo_level"] == level)]
        m = g.merge(c, on="geo_id", how="inner")
        for universe_name, sub in (("all_scored", m), ("rank_eligible", m[m["rank_eligible"].astype(bool)])):
            ok = sub["clinic_capacity"].notna() & sub["measurement_capacity"].notna()
            okp = sub["cc_impl_providers_in_geo_per_100k"].notna() & sub["active_clinicians_per_100k_adults"].notna()
            rho, p = spearmanr(sub.loc[ok, "clinic_capacity"], sub.loc[ok, "measurement_capacity"]) \
                if ok.sum() > 2 else (np.nan, np.nan)
            rho2, _ = spearmanr(sub.loc[okp, "cc_impl_providers_in_geo_per_100k"],
                                sub.loc[okp, "active_clinicians_per_100k_adults"]) if okp.sum() > 2 else (np.nan, None)
            top = sub[ok]
            k = max(1, int(round(0.1 * len(top))))
            ov = len(set(top.nlargest(k, "clinic_capacity")["geo_id"]) &
                     set(top.nlargest(k, "measurement_capacity")["geo_id"])) / k if len(top) else np.nan
            rows.append({"measurement_id": mid, "geo_level": level, "universe": universe_name,
                         "capacity_basis": c["capacity_basis"].iloc[0] if len(c) else "",
                         "n_geos": int(ok.sum()),
                         "spearman_measurement_capacity_vs_clinic_capacity": rho, "p_value": p,
                         "spearman_clinician_rate_vs_implementer_provider_rate": rho2,
                         "geos_with_zero_active_clinicians": int((sub["active_clinicians_n"].fillna(0) == 0).sum()),
                         "top_decile_n": k, "top_decile_overlap_share": ov})
    return pd.DataFrame(rows)


@lru_cache(maxsize=None)
def _clin_weights(measurement_id: str) -> np.ndarray | None:
    basis = capacity_basis(measurement_id)
    if basis == "none":
        return None
    w, _ = _npi_weights(provider_summary(measurement_id), basis)
    return w["clin"]


def slot_states() -> np.ndarray:
    """State FIPS of each slot NPI (for within-state permutations)."""
    slots, slot_codes, _ = location_slots()
    return slots["state_fips"].to_numpy()[slot_codes]


def capacity_pct_for(measurement_id: str, level: str, geo_ids, perm: np.ndarray | None = None) -> np.ndarray:
    """measurement_capacity percentile aligned to `geo_ids`, optionally with the slot NPIs' locations permuted
    (perm[j] = the slot NPI j is moved to). Used by the scoring (clinic_capacity_basis) and its Test 6b shuffle."""
    w = _clin_weights(measurement_id)
    if w is None:
        return np.full(len(geo_ids), np.nan)
    _, slot_codes, _ = location_slots()
    sc = slot_codes if perm is None else slot_codes[perm]
    comp = capacity_components(level, capacity_counts(level, {"clin": w}, sc))
    cap = pd.Series(comp["measurement_capacity"], index=universe(level)["geo_id"].to_numpy())
    return percentile_rank(cap.reindex(list(geo_ids)).to_numpy(float))


def location_shuffle(cap: pd.DataFrame, measurement_ids: list[str], draws: int = 200, state_draws: int = 100,
                     seed: int = SEED + 61) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Test-6b analogue: permute the location slots of all individual NPIs (nationally / within state), recompute
    county measurement_capacity, Spearman with the observed value. Activity stays with the NPI."""
    slots, slot_codes, npi_index = location_slots()
    st = slots["state_fips"].to_numpy()[slot_codes]
    by_state = [np.flatnonzero(st == s) for s in np.unique(st)]
    rng = np.random.default_rng(seed)
    perms = [("national", lambda: rng.permutation(len(slot_codes)), draws)]

    def within_state():
        p = np.arange(len(slot_codes))
        for idx in by_state:
            p[idx] = idx[rng.permutation(len(idx))]
        return p
    perms.append(("within_state", within_state, state_draws))
    weights = {}
    for mid in measurement_ids:
        basis = capacity_basis(mid)
        if basis == "none":
            continue
        summ = provider_summary(mid)
        w, _ = _npi_weights(summ, basis)
        weights[mid] = w["clin"]
    obs = {mid: capacity_components("county", capacity_counts("county", {"clin": w}, slot_codes))
           for mid, w in weights.items()}
    draws_rows = []
    for null, make, nd in perms:
        for d in range(nd):
            perm = make()
            sc = slot_codes[perm]  # NPI j occupies the slot of NPI perm[j]
            for mid, w in weights.items():
                comp = capacity_components("county", capacity_counts("county", {"clin": w}, sc))
                a, b = obs[mid]["measurement_capacity"], comp["measurement_capacity"]
                ok = ~np.isnan(a) & ~np.isnan(b)
                draws_rows.append({"measurement_id": mid, "null": null, "draw": d,
                                   "spearman_shuffled_vs_observed": spearmanr(a[ok], b[ok])[0]})
    dr = pd.DataFrame(draws_rows)
    summ = dr.groupby(["measurement_id", "null"])["spearman_shuffled_vs_observed"].agg(
        n_draws="size", mean="mean", p05=lambda s: s.quantile(0.05), p95=lambda s: s.quantile(0.95)).reset_index()
    summ["capacity_basis"] = summ["measurement_id"].map(capacity_basis)
    summ["reference_clinic_capacity_test6b_national"] = summ["measurement_id"].map(
        {"wearable_autonomic_activity_monitoring": 0.965, "nailfold_capillaroscopy": 0.953,
         "autonomic_function_testing": 0.937, "exercise_capacity_testing": 0.633})
    return summ, dr


# ------------------------------------------------------------------------------------------------------------------
# entry point
# ------------------------------------------------------------------------------------------------------------------

def report(draws: int = 200, state_draws: int = 100, shuffle: bool = True) -> dict:
    """Comparison with the scoring's clinic_capacity (needs deployment_opportunities) + the location shuffle.
    Reads geo_measurement_capacity; run it after the scoring step (pipeline step facilities_activity_report)."""
    cap = read_table(CAPACITY_TABLE)
    TABLES.mkdir(parents=True, exist_ok=True)
    cmp_ = compare_with_clinic_capacity(cap)
    cmp_.to_csv(TABLES / "stage_f_capacity_vs_clinic_capacity.csv", index=False)
    out = {"comparison": cmp_.to_dict("records")}
    if shuffle:
        bundles = [m for m, v in measurement_members().items() if v["kind"] == "bundle"]
        sh, dr = location_shuffle(cap, bundles, draws=draws, state_draws=state_draws)
        sh.to_csv(TABLES / "stage_f_capacity_location_shuffle.csv", index=False)
        dr.to_csv(TABLES / "stage_f_capacity_location_shuffle_draws.csv", index=False)
        out["shuffle"] = sh.to_dict("records")
    return out


def build() -> dict:
    """Tables only (no dependency on the scoring): geo_measurement_capacity, provider_measurement_summary,
    measurement_experience_sites. Pipeline step facilities_activity (before scoring)."""
    cap, summ = build_capacity()
    write_table(cap, CAPACITY_TABLE, producer=PRODUCER,
                description="Measurement-specific clinician capacity (Medicare FFS billing) per county/state")
    if len(summ):
        write_table(summ, SUMMARY_TABLE, producer=PRODUCER,
                    description="NPI x measurement id: Medicare FFS services by code role, codes billed, location")
    sites = build_experience_sites()
    write_table(sites, SITES_TABLE, producer=PRODUCER,
                description="U.S. facilities with registered trials / NIH projects mentioning a measurement class")
    return {"capacity_rows": int(len(cap)), "summary_rows": int(len(summ)), "experience_site_rows": int(len(sites))}


def run(draws: int = 200, state_draws: int = 100, shuffle: bool = True, with_report: bool = True) -> dict:
    t0 = time.time()
    out = build()
    if with_report:
        out.update(report(draws, state_draws, shuffle))
    out["seconds"] = round(time.time() - t0, 1)
    return out


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description="Measurement-specific capacity and experience sites (stage F)")
    ap.add_argument("--draws", type=int, default=200)
    ap.add_argument("--state-draws", type=int, default=100)
    ap.add_argument("--no-shuffle", action="store_true")
    ap.add_argument("--no-report", action="store_true", help="tables only (no comparison / shuffle)")
    a = ap.parse_args(argv)
    print(json.dumps(run(a.draws, a.state_draws, shuffle=not a.no_shuffle, with_report=not a.no_report), indent=2,
                     default=str))


if __name__ == "__main__":
    main()
