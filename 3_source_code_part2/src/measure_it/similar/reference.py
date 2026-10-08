"""Apply a computable phenotype to NHANES reference cohorts and estimate, design-based, how common the
device-positive-LIKE profile is among US adults (aggregate only).

What the numbers are, and are not
    NHANES 2011-2014 (`nhanes`) and 2003-2006 (`nhanes0306`) are pre-pandemic cross-sectional surveys of the US
    civilian non-institutionalised population. They contain no Long COVID, PTLDS or capillaroscopy item. The estimate
    is the survey-weighted share of adults whose OBSERVABLE features (demographics, self-reported history mapped to
    ICD-10-CM, labs mapped to LOINC, survey items crosswalked by stage B+C) put them on the positive side of the
    phenotype's decision rule. It is a "people who resemble the device-positive patients" profile prevalence in a
    general population, NOT a prevalence of Long COVID, of Lyme disease, or of the device finding.
    (The post-pandemic NHANES August 2021-August 2023 cycle asked post-COVID items, COQ.160-COQ.200, but the COQ_L data
    file is not released; see NHANES_2021_2023_FINDING.)

Method
    * Features: `measure_it.harmonize.features.feature_matrix([dataset])` (person_concepts); score:
      `measure_it.harmonize.phenotype.score_phenotype`. Features that do not exist in the reference, and item-missing
      values, are held at their centre (reduced phenotype; the OMOP and clinic packs use the same rule); coverage =
      share of |coefficient_standardized| mass that exists, overall and per domain; the weighted share of adults
      with every existing feature observed is reported. SURVEY features count as existing only through a
      `same_item` crosswalk (or the identical item).
    * Population: adults >= 18 examined in the MEC (MEC weights are valid for interview and examination items).
      When an existing laboratory feature comes from a random subsample (morning fasting: WTSAF2YR; 2011-2012
      subsample A: WTSA2YR) the subsample weight is read from the raw NHANES files and used instead (pooled over the
      cycles in which the variable exists). Weights, strata and PSUs are never invented: a dataset without design
      columns is skipped with a reason.
    * Estimator: `measure_it.similar.survey.domain_estimate` (Taylor linearisation, domain estimation over the full
      design; checked against `samplics` on NHANES: identical SE), 95% Wald-t CI on the design df and Korn-Graubard
      CI with the NCHS reliability rule.
    * NMR platform: a LOINC feature that the phenotype took from Nightingale NMR metabolomics is evaluated on NHANES
      routine-lab values; the two are not calibrated. When such features exist, a second row set
      (`variant="excluding_nmr_mapped"`) drops them.
    * Nearest-profile similarity: Gower distance of each reference adult to the device-positive centroid over the
      features shared by both, reported only as weighted quantiles.

Tables (aggregate, non-partition; rows for a phenotype are replaced on re-run):
    similar_reference_prevalence, similar_reference_profile, similar_reference_distance
"""
from __future__ import annotations

import json
import re
from dataclasses import dataclass

import numpy as np
import pandas as pd

from ..config import PROCESSED, UNKNOWN, raw_dir, utc_now_iso
from ..provenance import add_provenance
from ..store import write_table
from .phenotype_io import MIN_COVERAGE, coverage, coverage_flag, is_nmr, threshold_lp
from .survey import domain_estimate, weighted_quantiles

PRODUCER = "measure_it.similar.reference"
T_PREV, T_PROFILE, T_DIST = "similar_reference_prevalence", "similar_reference_profile", "similar_reference_distance"

NHANES_2021_2023_FINDING = (
    "NHANES August 2021-August 2023 (suffix _L) is published (DEMO_L etc., September 2024). Its CAPI COVID-19 "
    "questionnaire (COQ) asks post-COVID items (COQ.160 symptoms >= 4 weeks after infection, COQ.170 symptoms, "
    "COQ.180 activity limitation, COQ.190 still symptomatic, COQ.200 duration), but no COQ_L data file is released: "
    "it is absent from the Questionnaire data listing and the Limited Access (RDC) listing, and "
    "DataFiles/COQ_L.xpt returns HTTP 404 (checked 2026-10-07). RXQ_RX_L has no drug names or reason codes.")


@dataclass
class ReferenceSpec:
    dataset_id: str
    label: str
    participants_table: str
    labs_table: str
    raw_source: str
    cycles: dict                 # cycle_code -> {"subsample_files": {weight_var: file}}
    mec_weight: str = "wtmec4yr_pooled"
    age_topcode: int = 80


REFERENCES = {
    "nhanes": ReferenceSpec(
        "nhanes", "NHANES 2011-2014 (pooled 4-year MEC weights), US civilian non-institutionalised adults 18+; "
                  "pre-pandemic, no Long COVID item",
        "participants__nhanes", "participant_labs__nhanes", "nhanes_2011_2014",
        {"G": {"WTSAF2YR": "GLU_G", "WTSA2YR": "THYROD_G"}, "H": {"WTSAF2YR": "GLU_H"}}, age_topcode=80),
    "nhanes0306": ReferenceSpec(
        "nhanes0306", "NHANES 2003-2006 (pooled 4-year MEC weights), US civilian non-institutionalised adults 18+; "
                      "pre-pandemic, no Long COVID item",
        "participants__nhanes0306", "participant_labs__nhanes0306", "nhanes_2003_2006",
        {"C": {"WTSAF2YR": "L10AM_C"}, "D": {"WTSAF2YR": "GLU_D"}}, age_topcode=85),
}
DESIGN_COLS = ("sdmvstra", "sdmvpsu")
DIST_QS = [0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95]
MIN_CELL_N = 10   # unweighted n below which a stratum's distance quantiles are not reported


# --------------------------------------------------------------------------------------------- scoring


def transform_column(f: dict, x: pd.Series) -> pd.Series:
    """Transformed feature value with the stage-BC missing rule: a missing value sits at the feature centre
    (presence: the training prevalence `center`, default 0; standardized: contribution 0)."""
    t = f.get("transform", "presence")
    x = pd.to_numeric(x, errors="coerce").astype(float)
    if t == "presence":
        return x.fillna(float(f.get("center", 0.0)))
    if t == "standardized_value":
        return ((x - float(f["center"])) / (float(f["scale"]) or 1.0)).fillna(0.0)
    raise ValueError(f"unsupported transform {t!r} for {f['feature_key']}")


def linear_predictor(ph: dict, X: pd.DataFrame) -> pd.Series:
    """Local re-implementation of the stage-BC rule (tests only): intercept + sum(coef x t(x)); a feature absent from X
    sits at its centre for everyone."""
    lp = pd.Series(float(ph["intercept"]), index=X.index)
    for f in ph["features"]:
        x = X[f["feature_key"]] if f["feature_key"] in X.columns else pd.Series(np.nan, index=X.index)
        lp = lp + float(f["coefficient"]) * transform_column(f, x)
    return lp


def score(ph: dict, X: pd.DataFrame) -> pd.Series:
    """Phenotype-positive (bool) per row of X with `measure_it.harmonize.phenotype.score_phenotype` (stage BC).
    Features whose columns are not in X are at their centre (reduced phenotype)."""
    from ..harmonize.phenotype import score_phenotype
    return _positive_from_score(ph, score_phenotype(ph, X), X.index)


def _positive_from_score(ph: dict, s, index) -> pd.Series:
    thr = float(ph["decision_threshold"])
    if isinstance(s, pd.DataFrame):
        for c in ("positive", "phenotype_positive", "is_positive"):
            if c in s.columns:
                return s[c].astype(bool).reindex(index)
        for c in ("probability", "prob", "p", "score"):
            if c in s.columns:
                return (s[c] >= thr).reindex(index)
        for c in ("lp", "linear_predictor", "logit"):
            if c in s.columns:
                return (s[c] >= threshold_lp(ph)).reindex(index)
        raise ValueError(f"score_phenotype returned unrecognised columns {list(s.columns)}")
    s = pd.Series(s, index=index) if not isinstance(s, pd.Series) else s.reindex(index)
    if s.dtype == bool:
        return s
    return s >= thr


# --------------------------------------------------------------------------------------------- reference data


def load_design(spec: ReferenceSpec) -> pd.DataFrame:
    p = PROCESSED / f"{spec.participants_table}.parquet"
    if not p.exists():
        raise FileNotFoundError(f"{p} missing (run the NHANES ingestion)")
    d = pd.read_parquet(p)
    missing = [c for c in (*DESIGN_COLS, spec.mec_weight, "age_years", "sex", "mec_examined") if c not in d.columns]
    if missing:
        raise KeyError(f"{spec.participants_table} lacks design/demographic columns {missing}; not inventing them")
    return d[["participant_id", "seqn", "cycle_code", "age_years", "sex", "mec_examined", *DESIGN_COLS,
              spec.mec_weight]].copy()


def lab_subsample_weights(spec: ReferenceSpec) -> dict[str, str]:
    """NHANES lab variable -> subsample weight variable, parsed from the labs table's `subsample` text."""
    p = PROCESSED / f"{spec.labs_table}.parquet"
    if not p.exists():
        return {}
    l = pd.read_parquet(p, columns=["lab_variable", "subsample"]).drop_duplicates()
    out = {}
    for v, s in zip(l.lab_variable, l.subsample.astype(str)):
        m = re.search(r"\(weight (WT[A-Z0-9]+)\)", s)
        if m and m.group(1) not in ("WTMEC2YR", "WTINT2YR"):
            out[v] = m.group(1)
    return out


def read_subsample_weight(spec: ReferenceSpec, weight_var: str) -> pd.DataFrame:
    """seqn -> pooled subsample weight read from the raw NHANES files (divided by the number of cycles holding it)."""
    from ..ingestion import nhanes as N
    frames = []
    files = [(c, v[weight_var]) for c, v in spec.cycles.items() if weight_var in v]
    for cyc, f in files:
        path = raw_dir(spec.raw_source) / f"{f}.xpt"
        if not path.exists():
            raise FileNotFoundError(f"{path} missing; cannot read {weight_var}")
        df = N.read_xport(path, ["SEQN", weight_var])
        frames.append(df)
    w = pd.concat(frames, ignore_index=True)
    w["seqn"] = w["SEQN"].astype("int64")
    w["w"] = w[weight_var] / len(files)
    return w[["seqn", "w"]]


def feature_source_variables(dataset_id: str, ph_features: list[dict]) -> dict[str, set]:
    """feature_key -> NHANES source variables (from person_concepts__<dataset>.source_variable)."""
    p = PROCESSED / f"person_concepts__{dataset_id}.parquet"
    out: dict[str, set] = {}
    if not p.exists():
        return out
    pc = pd.read_parquet(p, columns=["vocabulary", "concept_code", "source_variable"]).drop_duplicates()
    for f in ph_features:
        if f["vocabulary"] != "LOINC":
            continue
        hit = pc[(pc.vocabulary == "LOINC") & pc.concept_code.isin(f["codes"])]
        out[f["feature_key"]] = {str(s).split(".")[-1].split(":")[-1] for s in hit.source_variable.dropna()}
    return out


def age_band(age: pd.Series, topcode: int) -> pd.Series:
    """Contract 10-year bands via the stage-BC helper (NHANES top-code decade -> '80+'); adults 18-19 -> '18-19'."""
    from ..harmonize.schema import age_band_from_years
    a = pd.to_numeric(pd.Series(age), errors="coerce")
    band = age_band_from_years(a, topcode)
    return band.where(~((a >= 18) & (a < 20)), "18-19")


# --------------------------------------------------------------------------------------------- analysis


def evaluable_features(ph: dict, X: pd.DataFrame) -> dict[str, bool]:
    return {f["feature_key"]: bool(f["feature_key"] in X.columns and X[f["feature_key"]].notna().any())
            for f in ph["features"]}


def device_positive_centroid(ph: dict) -> dict[str, float]:
    """Device-positive feature means carried by the phenotype JSON (aggregate), if stage B+C wrote them."""
    for k in ("device_positive_centroid", "device_positive_profile", "centroid"):
        c = ph.get(k)
        if isinstance(c, dict) and c:
            return {kk: float(v) for kk, v in c.items() if v is not None}
    prof = ph.get("feature_profile") or ph.get("feature_summaries")
    if isinstance(prof, list):
        out = {}
        for r in prof:
            m = r.get("mean_positive", r.get("device_positive_mean"))
            if r.get("feature_key") and m is not None:
                out[r["feature_key"]] = float(m)
        return out
    return {}


def gower_to_centroid(X: pd.DataFrame, centroid: dict[str, float], binary: set[str]) -> tuple[pd.Series, list[str]]:
    """Gower distance of each row to a centroid over shared features (missing values skipped per row).
    Numeric ranges are the reference's 1st-99th percentile range."""
    feats = [k for k in centroid if k in X.columns and X[k].notna().any()]
    num = pd.DataFrame(index=X.index)
    for k in feats:
        x = pd.to_numeric(X[k], errors="coerce")
        if k in binary:
            num[k] = (x.fillna(0) > 0).astype(float).sub(centroid[k]).abs()
        else:
            lo, hi = np.nanpercentile(x, [1, 99]) if x.notna().any() else (np.nan, np.nan)
            r = hi - lo
            num[k] = ((x.clip(lo, hi) - centroid[k]).abs() / r).clip(upper=1.0) if r and r > 0 else np.nan
    if not feats:
        return pd.Series(np.nan, index=X.index), feats
    return num.mean(axis=1, skipna=True), feats


def analyse_reference(ph: dict, X: pd.DataFrame, design: pd.DataFrame, spec: ReferenceSpec,
                      scorer=None, subsample_weight: pd.DataFrame | None = None,
                      subsample_label: str = "") -> dict[str, pd.DataFrame]:
    """Core estimation for one phenotype x one reference (pure: no I/O). X is indexed by participant_id."""
    scorer = scorer or score
    ev = evaluable_features(ph, X)
    d = design.set_index("participant_id")
    weight = d[spec.mec_weight].where(d["mec_examined"].astype(bool))
    wlabel = f"{spec.mec_weight} (MEC examined)"
    if subsample_weight is not None:
        weight = d["seqn"].map(subsample_weight.set_index("seqn")["w"]).astype(float)
        wlabel = subsample_label
    adult = (pd.to_numeric(d["age_years"], errors="coerce") >= 18)
    eligible = adult & weight.gt(0)
    Xa = X.reindex(d.index)
    band = age_band(d["age_years"], spec.age_topcode)
    variants = [("all_available_features", [f for f in ph["features"] if ev[f["feature_key"]]])]
    nmr = [f for f in ph["features"] if ev[f["feature_key"]] and is_nmr(f)]
    if nmr:
        variants.append(("excluding_nmr_mapped", [f for f in variants[0][1] if not is_nmr(f)]))
    # every feature not used in a variant is absent from the scored matrix -> at its centre (stage-BC rule)
    prev_rows, prof_rows, dist_rows = [], [], []
    for variant, feats in variants:
        ev_v = {f["feature_key"]: (f in feats) for f in ph["features"]}
        cov = coverage(ph, ev_v)
        cols = [f["feature_key"] for f in feats]
        pos = scorer(ph, Xa.loc[eligible, cols]).reindex(d.index).fillna(False).astype(bool)
        y = pos.astype(float).where(eligible)
        keys = [f["feature_key"] for f in feats]
        complete = Xa[keys].notna().all(axis=1).astype(float) if keys else pd.Series(1.0, index=d.index)
        complete_share = domain_estimate(complete.where(eligible), weight, d.sdmvstra, d.sdmvpsu, eligible,
                                         proportion=False).estimate
        strata = [("all", "all", eligible)]
        for b in sorted(band[eligible].dropna().unique()):
            for sx in ("female", "male"):
                strata.append((b, sx, eligible & band.eq(b) & d["sex"].eq(sx)))
        for sx in ("female", "male"):
            strata.append(("all", sx, eligible & d["sex"].eq(sx)))
        for b in sorted(band[eligible].dropna().unique()):
            strata.append((b, "all", eligible & band.eq(b)))
        cflag = coverage_flag(cov["coverage"])
        for ab, sx, dom in strata:
            e = domain_estimate(y.values, weight.values, d.sdmvstra.values, d.sdmvpsu.values, dom.values)
            if cflag.startswith("insufficient"):
                for a in ("estimate", "se", "ci_low", "ci_high", "kg_ci_low", "kg_ci_high"):
                    setattr(e, a, np.nan)
                e.nchs_reliability = "suppress: feature coverage < %.2f" % MIN_COVERAGE
            prev_rows.append({"coverage_flag": cflag,
                "variant": variant, "age_band": ab, "sex": sx, "prevalence": e.estimate, "se": e.se,
                "ci_low": e.ci_low, "ci_high": e.ci_high, "kg_ci_low": e.kg_ci_low, "kg_ci_high": e.kg_ci_high,
                "design_df": e.df, "n_unweighted": e.n_unweighted,
                "n_positive_unweighted": int((pos & dom).sum()) if not cflag.startswith("insufficient") else None,
                "weighted_population": e.sum_weights, "nchs_reliability": e.nchs_reliability,
                "feature_coverage": cov["coverage"], "n_features_used": len(feats), "n_features": cov["n_features"],
                "coverage_by_domain": json.dumps({k: round(v["evaluable_share_of_domain"], 4)
                                                  for k, v in cov["by_domain"].items()}),
                "weighted_share_complete_features": complete_share, "weight_variable": wlabel})
        if variant != "all_available_features":
            continue
        for f in ph["features"]:
            k = f["feature_key"]
            row = {"feature_key": k, "domain": f["domain"], "label": f.get("label", ""),
                   "coefficient": float(f["coefficient"]), "transform": f.get("transform", ""),
                   "unit": f.get("unit"), "available_in_reference": ev[k], "nmr_platform_mismatch": is_nmr(f) and ev[k],
                   "reference_crosswalk": "; ".join(f"{r.get('reference')} ({r.get('match_quality')})"
                                                    for r in f.get("reference_crosswalk") or [])}
            if ev[k]:
                x = pd.to_numeric(Xa[k], errors="coerce")
                if f.get("transform") == "presence":
                    x = (x.fillna(0) > 0).astype(float)
                for grp, dom in (("phenotype_positive", eligible & pos), ("others", eligible & ~pos)):
                    e = domain_estimate(x.values, weight.values, d.sdmvstra.values, d.sdmvpsu.values, dom.values,
                                        proportion=False)
                    row[f"mean_{grp}"], row[f"se_{grp}"], row[f"n_{grp}"] = e.estimate, e.se, e.n_unweighted
            prof_rows.append(row)
        cen = device_positive_centroid(ph)
        if cen:
            binary = {f["feature_key"] for f in ph["features"] if f.get("transform") == "presence"}
            dist, used = gower_to_centroid(Xa.loc[eligible], cen, binary)
            dist = dist.reindex(d.index)
            for grp, dom in (("all_adults", eligible), ("phenotype_positive", eligible & pos),
                             ("others", eligible & ~pos)):
                ok = dom & dist.notna()
                n = int(ok.sum())
                qs = weighted_quantiles(dist[ok], weight[ok], DIST_QS) if n >= MIN_CELL_N else [np.nan] * len(DIST_QS)
                m = domain_estimate(dist.values, weight.values, d.sdmvstra.values, d.sdmvpsu.values, ok.values,
                                    proportion=False) if n >= MIN_CELL_N else None
                dist_rows.append({"group": grp, "metric": "gower_to_device_positive_centroid",
                                  "n_features_shared": len(used), "features_shared": ",".join(used),
                                  "n_unweighted": n, "weighted_mean": m.estimate if m else np.nan,
                                  "weighted_mean_se": m.se if m else np.nan,
                                  **{f"q{int(q * 100):02d}": v for q, v in zip(DIST_QS, qs)}})
        else:
            dist_rows.append({"group": "all_adults", "metric": "gower_to_device_positive_centroid",
                              "n_features_shared": 0, "features_shared": "", "n_unweighted": 0,
                              "note": "no device-positive centroid in the computable phenotype JSON "
                                      "(key device_positive_centroid); distance not computed"})
    return {"prevalence": pd.DataFrame(prev_rows), "profile": pd.DataFrame(prof_rows),
            "distance": pd.DataFrame(dist_rows), "evaluable": ev}


def choose_weight(ph: dict, ev: dict, spec: ReferenceSpec) -> tuple[pd.DataFrame | None, str, list[str]]:
    """Subsample weight when an existing LOINC feature comes from a random NHANES subsample."""
    notes = []
    sub = lab_subsample_weights(spec)
    srcs = feature_source_variables(spec.dataset_id, [f for f in ph["features"] if ev.get(f["feature_key"])])
    needed = {}
    for k, vars_ in srcs.items():
        for v in vars_:
            if v in sub:
                needed.setdefault(sub[v], set()).add(k)
    if any(f["vocabulary"] == "LOINC" and ev.get(f["feature_key"]) for f in ph["features"]) and not srcs:
        notes.append("person_concepts source_variable not available: LOINC features assumed to be full-MEC labs")
    if not needed:
        return None, "", notes
    sizes = {}
    for wv in needed:
        try:
            w = read_subsample_weight(spec, wv)
            sizes[wv] = (int((w.w > 0).sum()), w)
        except FileNotFoundError as e:
            notes.append(str(e))
    if not sizes:
        return None, "", notes
    wv = min(sizes, key=lambda k: sizes[k][0])
    if len(needed) > 1:
        notes.append(f"features from several subsamples {sorted(needed)}; the smallest ({wv}) weight is used, an "
                     "approximation (NHANES has no weight for the intersection)")
    notes.append(f"subsample weight {wv} (pooled) used because of {sorted(needed[wv])}")
    return sizes[wv][1], f"{wv} pooled (subsample: {', '.join(sorted(needed[wv]))})", notes


def same_item_crosswalk(ph: dict) -> dict[str, str]:
    """SURVEY feature -> reference feature key when the crosswalk says `same_item` (values comparable as they are).
    same_construct / related_construct items are NOT substituted (different scales); they stay not evaluable."""
    out = {}
    for f in ph["features"]:
        for r in f.get("reference_crosswalk") or []:
            ref = r.get("reference")
            if r.get("match_quality") == "same_item" and ref:
                key = ref if ref.startswith("SURVEY:") else f"SURVEY:{ref}"
                if key != f["feature_key"]:
                    out[f["feature_key"]] = key
                break
    return out


def _feature_matrix(dataset_id: str, keys: list[str], xwalk: dict[str, str] | None = None) -> pd.DataFrame:
    from ..harmonize.features import feature_matrix
    xwalk = xwalk or {}
    want = sorted(set(keys) | set(xwalk.values()))
    X = feature_matrix([dataset_id], want)
    if "participant_id" in X.columns:
        X = X.set_index("participant_id")
    for k, ref in xwalk.items():
        if k not in X.columns or X[k].isna().all():
            if ref in X.columns:
                X[k] = X[ref]
    X = X.dropna(axis=1, how="all")
    return X[[k for k in keys if k in X.columns]]


def run_reference(ph: dict, dataset_id: str, *, log=print) -> dict:
    spec = REFERENCES[dataset_id]
    design = load_design(spec)
    xwalk = same_item_crosswalk(ph)
    X = _feature_matrix(dataset_id, [f["feature_key"] for f in ph["features"]], xwalk)
    ev = evaluable_features(ph, X)
    sw, slabel, notes = choose_weight(ph, ev, spec)
    res = analyse_reference(ph, X, design, spec, subsample_weight=sw, subsample_label=slabel)
    res["notes"] = notes + [f"{k} evaluated with the same-item reference {v}" for k, v in xwalk.items() if ev.get(k)]
    for k in ("prevalence", "profile", "distance"):
        df = res[k]
        df.insert(0, "reference_dataset", dataset_id)
        df.insert(0, "phenotype_id", ph["phenotype_id"])
        df["reference_population"] = spec.label
        df["condition_id"] = ph.get("condition_id", UNKNOWN)
        df["notes"] = "; ".join(notes)
    log(f"[similar] {ph['phenotype_id']} x {dataset_id}: coverage "
        f"{res['prevalence']['feature_coverage'].iloc[0]:.2f}, national prevalence "
        f"{res['prevalence']['prevalence'].iloc[0]:.4f}")
    return res


def object_ids(df: pd.DataFrame) -> pd.Series:
    """`<phenotype_id>|<reference>|<variant>|<age_band or feature or group>|<sex>` (traceable aggregate rows)."""
    parts = [df["phenotype_id"].astype(str), df["reference_dataset"].astype(str)]
    for c in ("variant", "age_band", "sex", "feature_key", "group"):
        if c in df.columns:
            parts.append(df[c].fillna("").astype(str))
    return "similar:" + pd.concat(parts, axis=1).agg("|".join, axis=1)


SIMILAR_TABLES = ("similar_reference_prevalence", "similar_reference_profile", "similar_reference_distance")


def drop_dataset(ds: str) -> list[str]:
    """Remove a user dataset's phenotype rows (`byod_<id>:...`) from the aggregate similar_reference_* tables; a table
    left with no rows is deleted. Called by `byod remove`."""
    from ..config import PROCESSED
    from ..store import read_table, table_exists, write_table
    out = []
    for t in SIMILAR_TABLES:
        if not table_exists(t):
            continue
        df = read_table(t)
        own = df["phenotype_id"].astype(str).str.startswith(f"{ds}:")
        if not own.any():
            continue
        if own.all():
            for suf in (".parquet", ".meta.json"):
                (PROCESSED / f"{t}{suf}").unlink(missing_ok=True)
        else:
            write_table(df[~own], t, producer="measure_it.byod.remove", description=f"{t} without {ds}")
        out.append(t)
    return out


def _provenance(df: pd.DataFrame, rid: str, notes: str) -> pd.DataFrame:
    return add_provenance(df, data_layer="derived", source_name="NHANES (CDC NCHS) x user computable phenotype",
                          source_version="participants__nhanes / __nhanes0306 + person_concepts partitions",
                          retrieved_at=utc_now_iso(), evidence_type="survey_estimate", source_record_id=rid,
                          source_geographic_resolution="national", provenance_notes=notes)


CAVEAT = ("Survey-weighted share of US adults (pre-pandemic NHANES) whose observable features fall on the positive "
          "side of a device-derived phenotype rule: a resemblance profile, NOT Long COVID / Lyme prevalence and not "
          "the device finding. Reduced phenotype when feature_coverage < 1. " + NHANES_2021_2023_FINDING)


def write_tables(results: list[dict], phenotype_ids: list[str]) -> None:
    """Upsert the aggregate rows of these phenotypes into the three non-partition tables."""
    for key, name, desc in (("prevalence", T_PREV, "Design-based prevalence of a device-positive-like profile in "
                                                     "NHANES reference cohorts (aggregate)"),
                            ("profile", T_PROFILE, "Feature means, phenotype-positive vs others, NHANES reference"),
                            ("distance", T_DIST, "Gower distance to the device-positive centroid: weighted quantiles")):
        new = pd.concat([r[key] for r in results], ignore_index=True) if results else pd.DataFrame()
        if new.empty:
            continue
        new["object_id"] = object_ids(new)
        new = _provenance(new, "object_id", CAVEAT)
        path = PROCESSED / f"{name}.parquet"
        if path.exists():
            old = pd.read_parquet(path)
            old = old[~old["phenotype_id"].isin(phenotype_ids)]
            new = pd.concat([old, new], ignore_index=True)
        write_table(new, name, producer=PRODUCER, description=desc)
