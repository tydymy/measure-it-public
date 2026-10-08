"""Measurement adapters for user-supplied device data.

TabularFeatureAdapter   the default for `device_<name>.csv|parquet` files that already hold features:
                        per-participant rows are used as they are; per-day rows (`day_index`) are summarised per
                        participant (mean, SD and number of days per feature). Nothing is fitted.
load_custom_adapter     a partner's own MeasurementAdapter (manifest `devices.<name>.adapter: {module, class}`), a
                        .py file inside the dataset folder, registered with measurements.adapters.register_adapter.

Both satisfy measure_it.measurements.adapter.MeasurementAdapter, so the downstream contract
(measurements.adapters.validate_description / measurement_context / scores_to_long) is the same as for the wearable
and capillaroscopy adapters. Scores are descriptive measurements of physiology, never diagnoses.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from ..config import load_config
from ..measurements.adapter import MeasurementAdapter, MeasurementDescription
from ..measurements.adapters import register_adapter, validate_description

ATYPICALITY = "tabular_feature_atypicality"


def _class_signals(measurement_class: str) -> list[str]:
    for m in load_config("measurements")["measurement_classes"]:
        if m["id"] == measurement_class:
            return [str(s) for s in m.get("signals", [])] or [measurement_class]
    return [measurement_class]


class TabularFeatureAdapter(MeasurementAdapter):
    """Pre-computed per-participant (or per-day) features -> one row per participant_id."""

    def __init__(self, device: str = "device", measurement_class: str = "accelerometry", dataset_id: str = "byod",
                 description: str = "") -> None:
        self.device = device
        self.measurement_class = measurement_class
        self.dataset_id = dataset_id
        self.description = description
        self._features: list[str] = []
        self._long = False

    def preprocess(self, raw_data) -> pd.DataFrame:
        df = raw_data.copy() if isinstance(raw_data, pd.DataFrame) else pd.read_csv(raw_data)
        if "participant_id" not in df.columns:
            raise ValueError("device data need a participant_id column")
        df["participant_id"] = df["participant_id"].astype(str)
        self._long = "day_index" in df.columns
        keys = ["participant_id"] + (["day_index"] if self._long else [])
        feats = [c for c in df.columns if c not in keys]
        for c in feats:
            df[c] = pd.to_numeric(df[c], errors="coerce").astype(float)
        if self._long:
            df["day_index"] = pd.to_numeric(df["day_index"], errors="coerce").astype(int)
        self._features = feats
        return df[keys + feats]

    def embed(self, raw_data) -> pd.DataFrame:
        df = self.preprocess(raw_data)
        if not self._long:
            return df.drop_duplicates("participant_id").reset_index(drop=True)
        g = df.groupby("participant_id", sort=True)
        mean = g[self._features].mean().add_suffix("__mean")
        sd = g[self._features].std(ddof=1).add_suffix("__sd")
        nd = g[self._features].count().add_suffix("__n_days")
        out = pd.concat([mean, sd, nd], axis=1)
        out = out[[c for f in self._features for c in (f"{f}__mean", f"{f}__sd", f"{f}__n_days")]]
        return out.reset_index()

    def phenotype_score(self, embedding: pd.DataFrame) -> pd.DataFrame:
        """Label-free atypicality: mean |within-dataset z-score| over the embedding features (per-day counts left
        out). A description of how unusual a participant's features are in this dataset, not a phenotype claim."""
        feats = [c for c in embedding.columns if c != "participant_id" and not c.endswith("__n_days")]
        X = embedding[feats].astype(float)
        sd = X.std(ddof=0).replace(0, np.nan)
        Z = (X - X.mean()) / sd
        out = pd.DataFrame({"participant_id": embedding["participant_id"].astype(str)})
        out[ATYPICALITY] = Z.abs().mean(axis=1, skipna=True).to_numpy()
        for c in feats:
            out[f"{ATYPICALITY}__{c}"] = Z[c].to_numpy()
        out[f"{ATYPICALITY}__n_components"] = Z.notna().sum(axis=1).to_numpy()
        return out

    def describe_measurement(self) -> MeasurementDescription:
        return MeasurementDescription(
            adapter_id=f"byod_tabular_{self.dataset_id}_{self.device}",
            modality=self.description or f"user-supplied device features ({self.device})",
            measurement_ids=[self.measurement_class], signals=_class_signals(self.measurement_class),
            input_unit_of_observation="participant x day" if self._long else "participant",
            embedding_features=list(self._features), phenotypes=[ATYPICALITY], public_training_data=[],
            requires_private_data=True, status="implemented",
            limitations=["User-supplied pre-computed features; the engine does not see raw signals.",
                         "tabular_feature_atypicality is a label-free within-dataset summary, not a phenotype or a "
                         "diagnosis; it has no phenotype_axes link (UNKNOWN by design)."])


def load_custom_adapter(root: Path, spec: dict, register_as: str) -> MeasurementAdapter:
    """Import a partner's adapter class from a .py file inside the dataset folder and register it."""
    path = (Path(root) / spec["module"]).resolve()
    if not str(path).startswith(str(Path(root).resolve())) or path.suffix != ".py" or not path.exists():
        raise ValueError(f"adapter module {spec['module']!r} must be a .py file inside the dataset folder")
    mod_name = f"measure_it_byod_adapter_{register_as}"
    s = importlib.util.spec_from_file_location(mod_name, path)
    mod = importlib.util.module_from_spec(s)
    sys.modules[mod_name] = mod
    s.loader.exec_module(mod)
    cls = getattr(mod, spec["class"], None)
    if cls is None or not (isinstance(cls, type) and issubclass(cls, MeasurementAdapter)):
        raise ValueError(f"{spec['class']!r} in {path.name} is not a measure_it MeasurementAdapter subclass")
    register_adapter(register_as, cls, overwrite=True)
    adapter = cls()
    return adapter


def adapter_for_device(root: Path, manifest: dict, device: str, dataset_id: str) -> tuple[MeasurementAdapter, str]:
    """(adapter instance, registry name) for one device file; the tabular default unless the manifest names one."""
    dev = (manifest.get("devices") or {}).get(device) or {}
    cls = dev.get("measurement_class") or manifest["measurement"]["class"]
    if dev.get("adapter"):
        name = dev["adapter"].get("register_as") or f"{dataset_id}_{device}"
        a = load_custom_adapter(root, dev["adapter"], name)
    else:
        name = f"{dataset_id}_{device}_tabular"
        desc = dev.get("description") or manifest["measurement"].get("device", "")
        register_adapter(name, lambda: TabularFeatureAdapter(device, cls, dataset_id, desc), overwrite=True)
        a = TabularFeatureAdapter(device, cls, dataset_id, desc)
    return a, name


def run_adapter(adapter: MeasurementAdapter, raw: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame,
                                                                        MeasurementDescription]:
    """(embedding, scores, description) with the contract checks the downstream layers rely on."""
    emb = adapter.embed(raw)
    if "participant_id" not in emb.columns or emb["participant_id"].duplicated().any():
        raise ValueError(f"{type(adapter).__name__}.embed must return one row per participant_id")
    feats = [c for c in emb.columns if c != "participant_id"]
    bad = [c for c in feats if not pd.api.types.is_numeric_dtype(emb[c])]
    if bad:
        raise ValueError(f"{type(adapter).__name__}.embed returned non-numeric columns {bad[:5]}")
    desc = adapter.describe_measurement()
    probs = validate_description(desc)
    if probs:
        raise ValueError(f"{type(adapter).__name__}.describe_measurement violates the contract: {probs}")
    scores = adapter.phenotype_score(emb)
    return emb, scores, desc
