"""Score a computable phenotype (docs/HARMONIZATION_CONTRACT.md section 2) on any feature matrix.

    fm = feature_matrix(["nhanes"], [f["feature_key"] for f in phenotype["features"]])
    p = score_phenotype(phenotype, fm)                 # P(device-positive-like) per participant
    p.attrs["coverage"]                                # per feature: share of participants with a value
    phenotype_coverage(phenotype, fm)                  # per participant: share of |coefficient| mass observed

Linear predictor:  intercept + sum over features of coefficient * t(x)
    transform 'presence'            t(x) = x (0/1)
    transform 'standardized_value'  t(x) = (x - center) / scale
Missing values are never filled silently. `missing="center"` (default) puts a missing feature at its center (the
case-population mean, i.e. a neutral contribution relative to the training population) and records, in
`attrs`, how often that happened per feature and how much coefficient mass was observed; `missing="nan"` returns NaN
for any participant missing any feature. Features absent from the matrix altogether (e.g. a SURVEY item that a
reference cohort never asked) count as missing for everyone and are reported in attrs["missing_features"].
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from ..config import PROJECT_ROOT

SCHEMA_PATH = PROJECT_ROOT / "templates" / "byod" / "computable_phenotype.schema.json"
TRANSFORMS = ("presence", "standardized_value")


def _contrib(f: dict, x: pd.Series) -> pd.Series:
    if f["transform"] == "presence":
        return float(f["coefficient"]) * x
    if f["transform"] == "standardized_value":
        scale = float(f["scale"]) or 1.0
        return float(f["coefficient"]) * (x - float(f["center"])) / scale
    raise ValueError(f"unknown transform {f['transform']!r} (one of {TRANSFORMS})")


def score_phenotype(phenotype: dict, feature_matrix: pd.DataFrame, *, missing: str = "center",
                    output: str = "probability") -> pd.Series:
    """P(subgroup) (or the linear predictor with output='logit') for every row of `feature_matrix`.

    Coverage is reported in `result.attrs`: 'coverage' {feature_key: share of rows observed}, 'missing_features'
    (keys absent from the matrix), 'coefficient_mass_available' (share of sum |coefficient_standardized| whose
    columns exist), 'missing_policy', 'n_rows_with_any_missing'.
    """
    if missing not in ("center", "nan"):
        raise ValueError("missing must be 'center' or 'nan'")
    feats = phenotype["features"]
    idx = feature_matrix.index
    lp = pd.Series(float(phenotype["intercept"]), index=idx, dtype=float)
    any_missing = pd.Series(False, index=idx)
    coverage, absent = {}, []
    w = np.array([abs(float(f.get("coefficient_standardized", f["coefficient"]))) for f in feats], float)
    avail = 0.0
    for f, wi in zip(feats, w):
        key = f["feature_key"]
        if key in feature_matrix.columns:
            x = pd.to_numeric(feature_matrix[key], errors="coerce").astype(float)
            avail += wi
        else:
            x = pd.Series(np.nan, index=idx)
            absent.append(key)
        miss = x.isna()
        coverage[key] = float(1 - miss.mean()) if len(x) else 0.0
        any_missing |= miss
        if miss.any():
            x = x.fillna(float(f.get("center", 0.0)))
        lp = lp + _contrib(f, x)
    if missing == "nan":
        lp = lp.where(~any_missing)
    out = lp if output == "logit" else 1.0 / (1.0 + np.exp(-lp))
    out.name = phenotype.get("phenotype_id", "phenotype")
    out.attrs = {"coverage": coverage, "missing_features": absent,
                 "coefficient_mass_available": float(avail / w.sum()) if w.sum() > 0 else 1.0,
                 "missing_policy": ("missing feature values set to the feature center (training-population mean); "
                                    "see coverage" if missing == "center" else "NaN for any row missing a feature"),
                 "n_rows_with_any_missing": int(any_missing.sum())}
    return out


def phenotype_coverage(phenotype: dict, feature_matrix: pd.DataFrame) -> pd.DataFrame:
    """Per row: n_features_observed, share of |coefficient_standardized| mass observed, and ehr_evaluable mass."""
    feats = phenotype["features"]
    w = pd.Series({f["feature_key"]: abs(float(f.get("coefficient_standardized", f["coefficient"]))) for f in feats})
    obs = pd.DataFrame({k: (feature_matrix[k].notna() if k in feature_matrix.columns
                            else pd.Series(False, index=feature_matrix.index)) for k in w.index})
    tot = float(w.sum()) or 1.0
    return pd.DataFrame({"n_features_observed": obs.sum(axis=1).astype(int),
                         "coefficient_mass_observed": obs.mul(w, axis=1).sum(axis=1) / tot},
                        index=feature_matrix.index)


def load_schema() -> dict:
    return json.loads(Path(SCHEMA_PATH).read_text())


def validate_phenotype(phenotype: dict) -> list[str]:
    """JSON-schema problems (templates/byod/computable_phenotype.schema.json); [] = valid."""
    import jsonschema
    v = jsonschema.Draft202012Validator(load_schema())
    return [f"{'/'.join(str(p) for p in e.path) or '(top level)'}: {e.message}" for e in v.iter_errors(phenotype)]
