"""Reading a computable phenotype (docs/HARMONIZATION_CONTRACT.md §2) and the pieces every stage-D consumer shares:
locating the JSON, feature domains, feature coverage, the decision threshold on the linear-predictor scale, and the
base-population codes for a condition.

Coverage = share of the phenotype's |coefficient_standardized| mass (|coefficient| when no standardized value is
given; the weighting of measure_it.harmonize.phenotype) whose features can be evaluated in a target (a reference
cohort, an OMOP CDM, a clinic EHR), overall and per domain. A feature that cannot be evaluated is held at its centre
(standardized: contribution 0; presence: coefficient x training prevalence), so the target applies a REDUCED
phenotype; coverage says how much of the model that is.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

from ..config import PROCESSED, RESULTS

RESULTS_ROOT = RESULTS / "byod"
PREFIX = "byod_"
DOMAINS = ("demographic", "condition", "measurement", "drug", "survey", "device", "omics")
VOCAB_DOMAIN = {"ICD10CM": "condition", "LOINC": "measurement", "RXNORM": "drug", "DEMOG": "demographic",
                "SURVEY": "survey", "DEVICE": "device", "OMICS": "omics"}
NMR_PLATFORMS = ("nmr_nightingale", "nmr")

# Base populations in scope (HARMONIZATION_CONTRACT §0). ICD-10-CM codes come from the phenotype JSON when given,
# otherwise from condition_registry; PTLDS has no ICD-10-CM code, so its base population is a Lyme-history proxy.
BASE_POPULATION_FALLBACK = {
    "long_covid": {"icd10cm": ["U09.9"], "label": "Long COVID (ICD-10-CM U09.9, in use since 2021-10-01)",
                   "proxy": False},
    "lyme_disease": {"icd10cm": ["A69.2"], "label": "Lyme disease (ICD-10-CM A69.2, A69.20-A69.29)", "proxy": False},
    "ptlds": {"icd10cm": ["A69.2"], "label": "PROXY: any Lyme disease code (A69.2x). Post-treatment Lyme disease "
                                            "syndrome has no ICD-10-CM code; this base population is people with a "
                                            "recorded Lyme diagnosis, most of whom do not have PTLDS", "proxy": True},
}


def manifest_id_of(dataset_id: str) -> str:
    return dataset_id[len(PREFIX):] if dataset_id.startswith(PREFIX) else dataset_id


def engine_id(dataset_id: str) -> str:
    return dataset_id if dataset_id.startswith(PREFIX) else PREFIX + dataset_id


def resolve_target(target: str | Path) -> tuple[str, Path]:
    """`<folder>` (with manifest.yaml) or `<dataset_id>` -> (manifest dataset id, results/byod/<id>/)."""
    p = Path(target)
    if p.is_dir() and (p / "manifest.yaml").exists():
        import yaml
        mid = yaml.safe_load((p / "manifest.yaml").read_text())["dataset_id"]
    else:
        mid = manifest_id_of(str(target))
    return mid, RESULTS_ROOT / mid


def phenotype_paths(results_dir: Path) -> list[Path]:
    """computable_phenotype.json, and any computable_phenotype*.json next to it (one per subgroup)."""
    return sorted(results_dir.glob("computable_phenotype*.json"))


def load_phenotype(path: Path) -> dict:
    ph = json.loads(Path(path).read_text())
    for k in ("phenotype_id", "intercept", "features", "decision_threshold"):
        if k not in ph:
            raise ValueError(f"{path}: computable phenotype lacks {k!r} (HARMONIZATION_CONTRACT §2)")
    for f in ph["features"]:
        f.setdefault("vocabulary", f["feature_key"].split(":", 1)[0])
        f.setdefault("domain", VOCAB_DOMAIN.get(f["vocabulary"], "other"))
        f.setdefault("codes", [f["feature_key"].split(":", 1)[1]])
        f.setdefault("transform", "presence")
        f.setdefault("label", f["feature_key"])
    return ph


MIN_COVERAGE = 0.25   # below: the reduced rule is mostly centre values; prevalence is suppressed (not interpretable)
LOW_COVERAGE = 0.50   # below: reported, flagged as a substantially reduced rule


def coverage_flag(cov: float) -> str:
    if not np.isfinite(cov) or cov < MIN_COVERAGE:
        return f"insufficient: coverage < {MIN_COVERAGE:.2f}, prevalence suppressed (rule dominated by centre values)"
    if cov < LOW_COVERAGE:
        return f"low: coverage < {LOW_COVERAGE:.2f}, a substantially reduced rule"
    return "adequate"


def coef_weight(f: dict) -> float:
    """|coefficient on the standardized scale| (as measure_it.harmonize.phenotype), else |coefficient|."""
    return abs(float(f.get("coefficient_standardized", f["coefficient"])))


def centre_term(f: dict) -> float:
    """Contribution of a feature held at its centre: presence -> coef x center (training prevalence, default 0);
    standardized -> 0. Matches measure_it.harmonize.phenotype.score_phenotype(missing='center')."""
    return float(f["coefficient"]) * float(f.get("center", 0.0)) if f.get("transform", "presence") == "presence" else 0.0


def reduced_intercept(ph: dict, evaluable: dict[str, bool]) -> float:
    """Intercept of the reduced phenotype: features the target cannot evaluate are held at their centre."""
    return float(ph["intercept"]) + sum(centre_term(f) for f in ph["features"] if not evaluable.get(f["feature_key"]))


def threshold_lp(ph: dict) -> float:
    """Decision threshold on the linear-predictor scale (logit of the probability threshold)."""
    t = float(ph["decision_threshold"])
    if not 0 < t < 1:
        raise ValueError("decision_threshold must be a probability in (0, 1)")
    return math.log(t / (1 - t))


def is_nmr(f: dict) -> bool:
    text = " ".join(str(f.get(k, "")) for k in ("platform", "notes", "mapping_note")).lower()
    return any(p in text for p in NMR_PLATFORMS) or "nightingale" in text


def coverage(ph: dict, evaluable: dict[str, bool]) -> dict:
    """Share of |coefficient| mass evaluable, overall and per domain. `evaluable[feature_key]` -> bool."""
    rows = []
    for f in ph["features"]:
        rows.append({"feature_key": f["feature_key"], "domain": f["domain"], "abs_coef": coef_weight(f),
                     "evaluable": bool(evaluable.get(f["feature_key"], False))})
    df = pd.DataFrame(rows, columns=["feature_key", "domain", "abs_coef", "evaluable"])
    total = df["abs_coef"].sum()
    out = {"coverage": float(df.loc[df.evaluable, "abs_coef"].sum() / total) if total > 0 else float("nan"),
           "n_features": int(len(df)), "n_evaluable": int(df["evaluable"].sum()), "by_domain": {}}
    for d, g in df.groupby("domain"):
        out["by_domain"][d] = {"n_features": int(len(g)), "n_evaluable": int(g["evaluable"].sum()),
                               "abs_coef_mass": float(g["abs_coef"].sum()),
                               "share_of_total_mass": float(g["abs_coef"].sum() / total) if total > 0 else float("nan"),
                               "evaluable_share_of_domain": float(g.loc[g.evaluable, "abs_coef"].sum() / g["abs_coef"].sum())
                               if g["abs_coef"].sum() > 0 else float("nan")}
    return out


def base_population(ph: dict) -> dict:
    """{'condition_id', 'icd10cm': [...], 'label', 'proxy'} for the phenotype's base population."""
    cid = ph.get("condition_id") or ""
    bp = ph.get("base_population") or {}
    codes = list(bp.get("icd10cm") or [])
    fb = BASE_POPULATION_FALLBACK.get(cid, {})
    proxy = bool(fb.get("proxy", False)) and not codes
    if not codes:
        reg = PROCESSED / "condition_registry.parquet"
        if reg.exists():
            r = pd.read_parquet(reg, columns=["canonical_condition_id", "icd10cm_codes"])
            hit = r.loc[r.canonical_condition_id == cid, "icd10cm_codes"]
            if len(hit) and hit.iloc[0] is not None and len(hit.iloc[0]):
                codes = [c for c in hit.iloc[0] if "." in c or len(c) == 3]
    if not codes:
        codes = list(fb.get("icd10cm", []))
    # keep the most specific codes once (U09 + U09.9 -> U09.9 is covered by the U09 prefix; keep both harmlessly)
    label = bp.get("description") or fb.get("label") or (f"{cid}: ICD-10-CM {', '.join(codes)}" if codes else cid)
    return {"condition_id": cid, "icd10cm": sorted(set(codes)), "label": label, "proxy": proxy}


def annotate_platform(ph: dict) -> dict:
    """Copy the `platform` of each LOINC feature from the user dataset's person_concepts partition when the phenotype
    JSON does not carry it (Nightingale NMR measures are flagged there by stage B+C)."""
    ds = str(ph.get("phenotype_id", "")).split(":", 1)[0]
    p = PROCESSED / f"person_concepts__{ds}.parquet"
    if not ds.startswith(PREFIX) or not p.exists():
        return ph
    try:
        pc = pd.read_parquet(p, columns=["vocabulary", "concept_code", "platform"])
    except Exception:   # older partition without the platform column
        return ph
    plat = pc[pc.vocabulary == "LOINC"].dropna(subset=["platform"]).groupby("concept_code")["platform"].first()
    for f in ph["features"]:
        if f.get("vocabulary") == "LOINC" and not f.get("platform"):
            hits = {plat[c] for c in f.get("codes", []) if c in plat.index}
            if hits:
                f["platform"] = sorted(hits)[0]
    return ph
