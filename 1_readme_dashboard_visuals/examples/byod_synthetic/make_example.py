"""Generate the SYNTHETIC bring-your-own-data tutorial dataset (seeded; nothing here is a real person).

    uv run python examples/byod_synthetic/make_example.py                # -> examples/byod_synthetic/synthetic_demo/
    uv run python examples/byod_synthetic/make_example.py --out DIR --n-cases 40 --n-controls 50 --fast

Every file, id and record is marked SYNTHETIC: participant ids are SYNTHETIC-0001..., the manifest sets
`synthetic: true` and `demo: true` (so the performance record never enters a default ranking; `measure-it byod
deploy --demo` includes it), titles and label definitions say SYNTHETIC, and README_SYNTHETIC.md explains the
simulation. The effect sizes are invented to make the walkthrough show a moderate, imperfect signal; they say nothing
about ME/CFS.

Story (simulated): ME/CFS cases vs healthy controls; a chest-strap ECG patch worn for up to 14 days (daily RMSSD,
SDNN, resting heart rate, upright minutes), a wrist actigraph (per-participant summaries), routine EHR labs (LOINC),
EHR diagnosis codes (ICD-10-CM; G93.32 is label-defining and is excluded from models by the engine) and a small
serum proteomics panel. The ECG patch goes through a partner-written MeasurementAdapter (adapters/), the actigraph
through the engine's TabularFeatureAdapter.

A second SYNTHETIC dataset mimics the layer structure of the MAESTRO study (make_maestro; `--kind maestro`):
cohorts Long COVID, persistent symptoms after Lyme (ptlds), acute Lyme and healthy; nailfold capillaroscopy, a ring
wearable (daily HR, HRV, SpO2, sleep, steps), PRO instruments (COMPASS-31, FSS-9, DSQ-PEM, PROMIS Cognitive,
FUNCAP27), a Nightingale-like NMR panel, an Olink-like protein panel and self-reported history; NO EHR files. A
microvascular latent factor drives the capillaroscopy score and, within the Long COVID cases, Raynaud history,
COMPASS-31 vasomotor score and GlycA: the planted signal that `measure-it byod subgroup` should recover. Every
number is invented; nothing describes MAESTRO's data or results.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

SEED = 20260928
HERE = Path(__file__).resolve().parent

ADAPTER = '''"""SYNTHETIC example of a partner-written MeasurementAdapter (daily ECG-patch summaries -> per-person features).

It shows the contract only: preprocess / embed / phenotype_score / describe_measurement. Registered by
`measure-it byod ingest` through the manifest (devices.synthetic_hrv_patch.adapter).
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from measure_it.measurements.adapter import MeasurementAdapter, MeasurementDescription


class SyntheticHrvPatchAdapter(MeasurementAdapter):
    FEATURES = ["log_rmssd_median", "rmssd_day_cv", "sdnn_median_ms", "resting_hr_median_bpm",
                "upright_minutes_median", "n_valid_days"]

    def preprocess(self, raw_data) -> pd.DataFrame:
        df = raw_data.copy()
        for c in ("rmssd_ms", "sdnn_ms", "resting_hr_bpm", "upright_minutes"):
            df[c] = pd.to_numeric(df[c], errors="coerce")
        # a day counts when RMSSD is physiologically plausible (5-250 ms) and resting HR is 30-150 bpm
        ok = df["rmssd_ms"].between(5, 250) & df["resting_hr_bpm"].between(30, 150)
        return df[ok]

    def embed(self, raw_data) -> pd.DataFrame:
        d = self.preprocess(raw_data)
        g = d.groupby("participant_id")
        out = pd.DataFrame({
            "log_rmssd_median": np.log(g["rmssd_ms"].median()),
            "rmssd_day_cv": g["rmssd_ms"].std(ddof=1) / g["rmssd_ms"].mean(),
            "sdnn_median_ms": g["sdnn_ms"].median(),
            "resting_hr_median_bpm": g["resting_hr_bpm"].median(),
            "upright_minutes_median": g["upright_minutes"].median(),
            "n_valid_days": g.size().astype(float)})
        return out.reset_index()

    def phenotype_score(self, embedding: pd.DataFrame) -> pd.DataFrame:
        x = embedding["log_rmssd_median"].astype(float)
        z = (x - x.mean()) / x.std(ddof=0)
        return pd.DataFrame({"participant_id": embedding["participant_id"],
                             "decreased_hrv": (-z).to_numpy(),          # a phenotype_axes axis id -> HPO:0031861
                             "decreased_hrv__log_rmssd_z": z.to_numpy(),
                             "decreased_hrv__n_components": 1})

    def describe_measurement(self) -> MeasurementDescription:
        return MeasurementDescription(
            adapter_id="synthetic_hrv_patch_v1", modality="SYNTHETIC chest-strap ECG patch (daily HRV summaries)",
            measurement_ids=["hrv"], signals=["RMSSD", "SDNN", "resting heart rate", "upright time"],
            input_unit_of_observation="participant x day", embedding_features=list(self.FEATURES),
            phenotypes=["decreased_hrv"], public_training_data=[], requires_private_data=True,
            status="implemented",
            limitations=["SYNTHETIC example adapter; the within-dataset reference is the whole dataset (label-free)."])
'''


def make(out: Path, n_cases: int = 90, n_controls: int = 110, fast: bool = False, seed: int = SEED,
         dataset_id: str = "synthetic_demo") -> Path:
    rng = np.random.default_rng(seed)
    out.mkdir(parents=True, exist_ok=True)
    (out / "adapters").mkdir(exist_ok=True)
    n = n_cases + n_controls
    ids = [f"SYNTHETIC-{i:04d}" for i in range(1, n + 1)]
    y = np.r_[np.ones(n_cases, int), np.zeros(n_controls, int)]
    perm = rng.permutation(n)
    y = y[perm]
    female = rng.random(n) < np.where(y == 1, 0.75, 0.6)
    age = np.clip(rng.normal(np.where(y == 1, 43, 41), 11), 20, 79.9)
    band = [f"{int(a // 10) * 10}-{int(a // 10) * 10 + 9}" for a in age]
    pd.DataFrame({"participant_id": ids, "me_cfs": y, "age_band": band,
                  "sex": np.where(female, "female", "male")}).to_csv(out / "participants.csv", index=False)
    # latent person-level autonomic signal (shared by the ECG patch and, weakly, the actigraph)
    lat = rng.normal(0, 1, n) + np.where(y == 1, 1.0, 0.0)
    age_c = (age - 42) / 11
    rows = []
    for i, pid in enumerate(ids):
        days = [d for d in range(14) if rng.random() > 0.12]
        base_rmssd = np.exp(3.6 - 0.30 * lat[i] - 0.15 * age_c[i] + rng.normal(0, 0.18))
        base_hr = 64 + 3.5 * lat[i] + rng.normal(0, 4)
        base_up = 420 - 45 * lat[i] + rng.normal(0, 70)
        for d in days:
            rows.append({"participant_id": pid, "day_index": d,
                         "rmssd_ms": round(float(base_rmssd * np.exp(rng.normal(0, 0.22))), 1),
                         "sdnn_ms": round(float(base_rmssd * 1.35 * np.exp(rng.normal(0, 0.18))), 1),
                         "resting_hr_bpm": round(float(base_hr + rng.normal(0, 3)), 1),
                         "upright_minutes": round(float(max(20, base_up + rng.normal(0, 60))), 0)})
    pd.DataFrame(rows).to_csv(out / "device_synthetic_hrv_patch.csv", index=False)
    act = pd.DataFrame({"participant_id": ids,
                        "steps_mean_daily": np.round(np.exp(8.7 - 0.12 * lat + rng.normal(0, 0.35, n))),
                        "sedentary_fraction": np.round(np.clip(0.62 + 0.03 * lat + rng.normal(0, 0.07, n), 0.2, 0.95), 3),
                        "interdaily_stability": np.round(np.clip(0.55 - 0.02 * lat + rng.normal(0, 0.1, n), 0.05, 0.95), 3)})
    act = act.drop(index=rng.choice(n, size=int(0.08 * n), replace=False))   # 8% without an actigraph
    act.to_csv(out / "device_synthetic_wrist_actigraphy.csv", index=False)
    labs = []
    spec = [("1988-5", "C reactive protein", "mg/L", lambda k: np.exp(rng.normal(0.3 + 0.25 * y[k], 0.9))),
            ("718-7", "Hemoglobin", "g/dL", lambda k: rng.normal(13.6 - 0.9 * female[k], 1.1)),
            ("2276-4", "Ferritin", "ng/mL", lambda k: np.exp(rng.normal(4.3 - 0.5 * female[k], 0.7))),
            ("3016-3", "Thyrotropin", "m[IU]/L", lambda k: np.exp(rng.normal(0.6, 0.5))),
            ("2345-7", "Glucose", "mg/dL", lambda k: rng.normal(92, 11))]
    for k, pid in enumerate(ids):
        for code, name, unit, draw in spec:
            if rng.random() < 0.1:
                continue
            for d in range(1 + int(rng.random() < 0.3)):
                labs.append({"participant_id": pid, "loinc": code, "value": round(float(draw(k)), 2), "unit": unit,
                             "lab_name": name, "day_index": int(rng.integers(0, 30))})
    pd.DataFrame(labs).to_csv(out / "ehr_labs.csv", index=False)
    dx = []
    codes = [("G93.32", 0.70, 0.0), ("R53.83", 0.55, 0.06), ("F41.9", 0.20, 0.10), ("K58.9", 0.15, 0.08),
             ("E03.9", 0.08, 0.08), ("I10", 0.14, 0.15), ("G90.A", 0.10, 0.01), ("J45.909", 0.09, 0.08)]
    for k, pid in enumerate(ids):
        for code, p1, p0 in codes:
            if rng.random() < (p1 if y[k] else p0):
                dx.append({"participant_id": pid, "icd10cm": code})
    pd.DataFrame(dx).to_csv(out / "ehr_diagnoses.csv", index=False)
    # EHR medication list (RxNorm ingredient CUIs, checked against RxNav on 2026-10-07)
    meds = [("6918", 0.10, 0.03), ("8787", 0.12, 0.02), ("4452", 0.06, 0.0), ("6963", 0.07, 0.0),
            ("36437", 0.18, 0.08), ("10582", 0.08, 0.07), ("5640", 0.30, 0.15), ("6809", 0.04, 0.05)]
    rx = [{"participant_id": pid, "rxnorm": cui} for k, pid in enumerate(ids) for cui, p1, p0 in meds
          if rng.random() < (p1 if y[k] else p0)]
    pd.DataFrame(rx).to_csv(out / "ehr_medications.csv", index=False)
    prot = {"participant_id": ids}
    shift = {3: 0.45, 11: 0.35, 27: -0.4, 42: 0.3, 60: -0.3}
    for j in range(80):
        prot[f"PROT{j + 1:03d}"] = np.round(rng.normal(5 + 0.02 * j + shift.get(j, 0) * y, 1.0), 4)
    om = pd.DataFrame(prot)
    om = om.drop(index=rng.choice(n, size=int(0.3 * n), replace=False))      # proteomics on 70% of participants
    om.to_csv(out / "omics_synthetic_proteomics.csv", index=False)
    (out / "adapters" / "hrv_patch_adapter.py").write_text(ADAPTER)
    manifest = {
        "schema_version": 1, "dataset_id": dataset_id,
        "title": "SYNTHETIC tutorial dataset: ME/CFS vs healthy, ECG patch + actigraph + EHR + proteomics",
        "owner": "measure-it examples (SYNTHETIC)",
        "description": "SYNTHETIC data generated by examples/byod_synthetic/make_example.py (seed %d); no real people." % seed,
        "licence": "CC0 1.0 (SYNTHETIC data generated by examples/byod_synthetic/make_example.py)",
        "synthetic": True, "demo": True,
        "data_use": {"deidentified": True, "use_permitted": True,
                     "statement": "SYNTHETIC: generated data, no real person; no IRB or data-use agreement applies.",
                     "irb_or_dua_reference": "not applicable (SYNTHETIC)"},
        "labels": [{"column": "me_cfs", "condition": "ME/CFS",
                    "definition": "SYNTHETIC label: simulated ME/CFS case status (Canadian Consensus Criteria in the "
                                  "simulated story)", "label_basis": "clinical_case_definition"}],
        "comparator": {"type": "healthy", "description": "SYNTHETIC healthy controls"},
        "measurement": {"class": "hrv", "bundle": "autonomic_function_testing",
                        "device": "SYNTHETIC chest-strap ECG patch (daily HRV, resting heart rate, upright time)"},
        "devices": {"synthetic_hrv_patch": {"measurement_class": "hrv",
                                            "description": "SYNTHETIC chest-strap ECG patch",
                                            "adapter": {"module": "adapters/hrv_patch_adapter.py",
                                                        "class": "SyntheticHrvPatchAdapter",
                                                        "register_as": "synthetic_hrv_patch"}},
                    "synthetic_wrist_actigraphy": {"measurement_class": "accelerometry",
                                                   "description": "SYNTHETIC wrist actigraph summaries"}},
        "id_hashing": "none",
        "primary_analysis": {"label": "me_cfs", "feature_blocks": ["device_synthetic_hrv_patch"],
                             "combined_blocks": ["device_synthetic_hrv_patch", "device_synthetic_wrist_actigraphy",
                                                 "ehr_labs"],
                             "specificity": 0.90, "covariates": ["age_band", "sex"], "n_splits": 5,
                             "n_repeats": 3 if fast else 10, "n_permutations": 20 if fast else 200,
                             "n_bootstrap": 200 if fast else 1000},
    }
    (out / "manifest.yaml").write_text("# SYNTHETIC tutorial manifest (generated; edit the generator, not this file)\n"
                                       + yaml.safe_dump(manifest, sort_keys=False, width=110))
    (out / "README_SYNTHETIC.md").write_text(
        "# SYNTHETIC tutorial dataset\n\nGenerated by `examples/byod_synthetic/make_example.py` (seed "
        f"{seed}). Nothing here describes a real person. The simulated effect sizes were chosen so that the "
        "walkthrough (docs/WALKTHROUGH.md) shows a moderate, imperfect signal; they carry no information about "
        "ME/CFS. The manifest marks the dataset `synthetic: true` and `demo: true`: its performance record is left "
        "out of every ranking unless `measure-it byod deploy --demo` is run, and `measure-it byod remove "
        "synthetic_demo` removes everything it created.\n")
    return out


MAESTRO_COHORTS = {"long_covid": 110, "ptlds": 60, "lyme_disease": 40, "healthy": 70}


def make_maestro(out: Path, cohorts: dict | None = None, fast: bool = False, seed: int = SEED,
                 dataset_id: str = "synthetic_maestro") -> Path:
    """SYNTHETIC dataset with MAESTRO-like layers and a planted device-defined subgroup (see the module docstring)."""
    rng = np.random.default_rng(seed + 1)
    cohorts = dict(cohorts or MAESTRO_COHORTS)
    out.mkdir(parents=True, exist_ok=True)
    coh = np.concatenate([[c] * n for c, n in cohorts.items()])
    coh = coh[rng.permutation(len(coh))]
    n = len(coh)
    ids = [f"SYNTHETIC-M{i:04d}" for i in range(1, n + 1)]
    sick = np.isin(coh, ["long_covid", "ptlds"]).astype(float) + 0.5 * (coh == "lyme_disease")
    female = rng.random(n) < np.where(coh == "healthy", 0.55, 0.68)
    age = np.clip(rng.normal(np.where(coh == "lyme_disease", 47, 42), 12), 19, 84.9)
    band = [f"{int(a // 10) * 10}-{int(a // 10) * 10 + 9}" for a in age]
    lab = {c: np.where(coh == c, "1", np.where(coh == "healthy", "0", "")) for c in ("long_covid", "ptlds",
                                                                                     "lyme_disease")}
    pd.DataFrame({"participant_id": ids, **lab, "age_band": band,
                  "sex": np.where(female, "female", "male")}).to_csv(out / "participants.csv", index=False)
    # latent microvascular factor (capillaroscopy) and disease burden
    m = rng.normal(0, 1, n) + np.select([coh == "long_covid", coh == "ptlds", coh == "healthy"], [0.3, 0.1, -0.8], 0.0)
    has_cap = rng.random(n) < 0.85
    cap = pd.DataFrame({
        "participant_id": ids,
        "abnormality_score": np.round(np.clip(1.5 + 0.9 * m + rng.normal(0, 0.45, n), 0, 4), 2),
        "capillary_density_per_mm": np.round(8.5 - 0.7 * m + rng.normal(0, 0.8, n), 2),
        "giant_capillaries_n": rng.poisson(np.exp(-1.0 + 0.6 * m)).astype(float),
        "microhemorrhages_n": rng.poisson(np.exp(-1.2 + 0.5 * m)).astype(float),
        "tortuosity_score": np.round(np.clip(1.2 + 0.5 * m + rng.normal(0, 0.6, n), 0, 3), 2)})[has_cap]
    cap.to_csv(out / "device_capillaroscopy.csv", index=False)
    ring = []
    for i, pid in enumerate(ids):
        if rng.random() < 0.1:
            continue
        hr0 = 63 + 6 * sick[i] + rng.normal(0, 5)
        hrv0 = np.exp(3.7 - 0.25 * sick[i] + rng.normal(0, 0.3))
        for d in range(14):
            if rng.random() < 0.15:
                continue
            ring.append({"participant_id": pid, "day_index": d,
                         "hr_mean_bpm": round(float(hr0 + rng.normal(0, 3)), 1),
                         "hrv_rmssd_ms": round(float(hrv0 * np.exp(rng.normal(0, 0.2))), 1),
                         "spo2_mean_pct": round(float(np.clip(96.5 + rng.normal(0, 0.8), 88, 100)), 1),
                         "sleep_minutes": round(float(max(120, 420 + 20 * sick[i] + rng.normal(0, 50))), 0),
                         "steps": round(float(max(200, np.exp(8.9 - 0.45 * sick[i] + rng.normal(0, 0.4))))), })
    pd.DataFrame(ring).to_csv(out / "device_ring.csv", index=False)
    lc = (coh == "long_covid").astype(float)
    # PRO instruments (COMPASS-31 domains weighted as in the instrument's 0-100 total; values invented)
    vaso = np.clip(np.round(0.8 + 0.9 * sick + 0.55 * lc * m + rng.normal(0, 0.6, n), 2), 0, 4.17)
    orth = np.clip(np.round(14 * sick + rng.normal(4, 6, n), 1), 0, 40)
    secr = np.clip(np.round(4 * sick + rng.normal(2, 3, n), 1), 0, 12.86)
    gi = np.clip(np.round(8 * sick + rng.normal(3, 4, n), 1), 0, 21.43)
    bladder = np.clip(np.round(2 * sick + rng.normal(1, 1.5, n), 1), 0, 11.11)
    pupil = np.clip(np.round(2 * sick + rng.normal(1, 1, n), 1), 0, 10.0)
    comp = pd.DataFrame({"participant_id": ids, "orthostatic": orth, "vasomotor": vaso, "secretomotor": secr,
                         "gastrointestinal": gi, "bladder": bladder, "pupillomotor": pupil,
                         "total": np.round(orth + vaso + secr + gi + bladder + pupil, 1)})
    comp[rng.random(n) > 0.05].to_csv(out / "survey_compass31.csv", index=False)
    pd.DataFrame({"participant_id": ids, "total": np.round(np.clip(2.5 + 2.6 * sick + rng.normal(0, 1, n), 1, 7), 2)}
                 )[rng.random(n) > 0.05].to_csv(out / "survey_fss9.csv", index=False)
    pd.DataFrame({"participant_id": ids,
                  "frequency_score": np.round(np.clip(15 + 45 * sick + rng.normal(0, 15, n), 0, 100), 1),
                  "severity_score": np.round(np.clip(12 + 40 * sick + rng.normal(0, 15, n), 0, 100), 1)}
                 )[rng.random(n) > 0.05].to_csv(out / "survey_dsqpem.csv", index=False)
    pd.DataFrame({"participant_id": ids, "t_score": np.round(52 - 9 * sick + rng.normal(0, 6, n), 1)}
                 )[rng.random(n) > 0.05].to_csv(out / "survey_promis_cog.csv", index=False)
    pd.DataFrame({"participant_id": ids, "total": np.round(np.clip(5.5 - 2.2 * sick + rng.normal(0, 0.8, n), 0, 6), 2)}
                 )[rng.random(n) > 0.05].to_csv(out / "survey_funcap27.csv", index=False)
    # Nightingale-like NMR panel (Nightingale reporting units: mmol/L, g/L, umol/L)
    ng = pd.DataFrame({
        "participant_id": ids,
        "Total_C": np.round(rng.normal(5.0, 0.9, n), 3), "HDL_C": np.round(rng.normal(1.5, 0.35, n), 3),
        "Clinical_LDL_C": np.round(rng.normal(3.0, 0.75, n), 3), "LDL_C": np.round(rng.normal(1.9, 0.5, n), 3),
        "Total_TG": np.round(np.exp(rng.normal(0.2, 0.4, n)), 3), "Glucose": np.round(rng.normal(5.2, 0.5, n), 3),
        "Creatinine": np.round(rng.normal(75, 14, n), 2), "Albumin": np.round(rng.normal(42, 2.8, n), 2),
        "GlycA": np.round(0.78 + 0.03 * sick + 0.07 * lc * m + rng.normal(0, 0.06, n), 4),
        "ApoB": np.round(rng.normal(0.95, 0.2, n), 3), "ApoA1": np.round(rng.normal(1.55, 0.22, n), 3),
        "Lactate": np.round(np.exp(rng.normal(0.3, 0.3, n)), 3), "Citrate": np.round(rng.normal(0.07, 0.012, n), 4),
        "Ile": np.round(rng.normal(0.055, 0.012, n), 4), "Leu": np.round(rng.normal(0.09, 0.018, n), 4),
        "Val": np.round(rng.normal(0.22, 0.04, n), 4), "XL_VLDL_P": np.round(np.exp(rng.normal(-1, 0.6, n)), 4)})
    ng[rng.random(n) < 0.9].to_csv(out / "omics_nightingale.csv", index=False)
    ol = {"participant_id": ids}
    shift = {4: 0.4, 17: 0.35, 33: -0.3}
    for j in range(60):
        ol[f"OLINK{j + 1:03d}"] = np.round(rng.normal(4 + 0.03 * j + shift.get(j, 0.0) * sick, 1.0), 4)
    pd.DataFrame(ol)[rng.random(n) < 0.8].to_csv(out / "omics_olink.csv", index=False)
    # self-reported medical history (free-text names as a study form would hold them; some intentionally unlisted)
    hist = []
    for i, pid in enumerate(ids):
        def add(text, p):
            if rng.random() < p:
                hist.append({"participant_id": pid, "condition": text})
        add("Raynaud's phenomenon", 1 / (1 + np.exp(-(-2.6 + 1.6 * m[i] * lc[i] + 0.6 * sick[i]))))
        add("Long COVID", 0.75 * lc[i])
        add("COVID-19", 0.9 if lc[i] else 0.45)
        add("Lyme disease", 0.9 if coh[i] in ("ptlds", "lyme_disease") else 0.03)
        add("POTS", 0.05 + 0.2 * sick[i])
        add("migraine", 0.12 + 0.15 * sick[i])
        add("hypothyroidism", 0.07)
        add("anxiety", 0.12 + 0.08 * sick[i])
        add("asthma", 0.09)
        add("seasonal allergies", 0.25)                    # not in the curated map: stays unmapped (counted)
    pd.DataFrame(hist).to_csv(out / "self_report_history.csv", index=False)
    manifest = {
        "schema_version": 1, "dataset_id": dataset_id,
        "title": "SYNTHETIC tutorial dataset: MAESTRO-like layers (capillaroscopy, ring, PROs, NMR, Olink, history)",
        "owner": "measure-it examples (SYNTHETIC)",
        "description": ("SYNTHETIC data generated by examples/byod_synthetic/make_example.py --kind maestro (seed %d); "
                        "mimics the layer structure of an infection-associated chronic illness study; no real people."
                        % seed),
        "licence": "CC0 1.0 (SYNTHETIC data generated by examples/byod_synthetic/make_example.py)",
        "synthetic": True, "demo": True,
        "data_use": {"deidentified": True, "use_permitted": True,
                     "statement": "SYNTHETIC: generated data, no real person; no IRB or data-use agreement applies.",
                     "irb_or_dua_reference": "not applicable (SYNTHETIC)"},
        "labels": [
            {"column": "long_covid", "condition": "long_covid", "label_basis": "clinical_diagnosis",
             "definition": "SYNTHETIC label: simulated Long COVID cohort vs simulated healthy controls"},
            {"column": "ptlds", "condition": "ptlds", "label_basis": "clinical_case_definition",
             "definition": "SYNTHETIC label: simulated persistent symptoms after treated Lyme vs healthy controls"},
            {"column": "lyme_disease", "condition": "lyme_disease", "label_basis": "clinical_diagnosis",
             "definition": "SYNTHETIC label: simulated acute Lyme cohort vs healthy controls"}],
        "comparator": {"type": "healthy", "description": "SYNTHETIC healthy controls"},
        "measurement": {"class": "capillaroscopy", "device": "SYNTHETIC nailfold capillaroscopy (abnormality score)"},
        "devices": {"capillaroscopy": {"measurement_class": "capillaroscopy",
                                       "description": "SYNTHETIC nailfold capillaroscopy summary features"},
                    "ring": {"measurement_class": "wearable_heart_rate",
                             "description": "SYNTHETIC ring wearable, daily HR / HRV / SpO2 / sleep / steps"}},
        "omics": {"nightingale": {"platform": "nmr_nightingale", "units": "nightingale_standard",
                                  "description": "SYNTHETIC Nightingale-like blood NMR panel (Nightingale units)"},
                  "olink": {"platform": "olink", "description": "SYNTHETIC Olink-like NPX panel"}},
        "id_hashing": "none",
        "primary_analysis": {"label": "long_covid", "feature_blocks": ["device_capillaroscopy"],
                             "combined_blocks": ["device_capillaroscopy", "device_ring"],
                             "specificity": 0.90, "covariates": ["age_band", "sex"], "n_splits": 5,
                             "n_repeats": 3 if fast else 10, "n_permutations": 20 if fast else 200,
                             "n_bootstrap": 200 if fast else 1000},
        "subgroup_analysis": {
            "name": "capillaroscopy_abnormal",
            "description": "SYNTHETIC: Long COVID cases with an abnormal nailfold capillaroscopy score",
            "within_label": "long_covid",
            "device_positive": {"block": "device_capillaroscopy", "feature": "abnormality_score", "threshold": 2.0,
                                "direction": ">="},
            "ehr_domains": ["condition", "measurement", "survey", "demographic"],
            "max_features": 8, "n_splits": 5, "n_repeats": 3 if fast else 10,
            "n_permutations": 20 if fast else 200, "n_bootstrap": 200 if fast else 1000},
    }
    (out / "manifest.yaml").write_text("# SYNTHETIC tutorial manifest (generated; edit the generator, not this file)\n"
                                       + yaml.safe_dump(manifest, sort_keys=False, width=110))
    (out / "README_SYNTHETIC.md").write_text(
        "# SYNTHETIC MAESTRO-like tutorial dataset\n\nGenerated by `examples/byod_synthetic/make_example.py --kind "
        f"maestro` (seed {seed}). Nothing here describes a real person or the MAESTRO study's data: only its layer "
        "structure is mimicked (no EHR extract; PROs, self-reported history and an NMR panel are the bridge to "
        "public data). A latent microvascular factor drives the capillaroscopy abnormality score and, inside the "
        "Long COVID cases, Raynaud history, the COMPASS-31 vasomotor score and GlycA; `measure-it byod subgroup` "
        "should find that planted signal. `synthetic: true`, `demo: true`.\n")
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", type=Path, default=HERE / "synthetic_demo")
    ap.add_argument("--n-cases", type=int, default=90)
    ap.add_argument("--n-controls", type=int, default=110)
    ap.add_argument("--fast", action="store_true", help="fewer CV repeats, permutations and bootstrap resamples")
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--dataset-id", default=None)
    ap.add_argument("--kind", choices=["ehr", "maestro"], default="ehr",
                    help="ehr: ME/CFS + ECG patch + EHR files (default); maestro: MAESTRO-like layers, no EHR")
    a = ap.parse_args()
    if a.kind == "maestro":
        out = make_maestro(a.out if a.out != HERE / "synthetic_demo" else HERE / "synthetic_maestro", fast=a.fast,
                           seed=a.seed, dataset_id=a.dataset_id or "synthetic_maestro")
        print(f"SYNTHETIC MAESTRO-like dataset written to {out}")
        return
    out = make(a.out, a.n_cases, a.n_controls, a.fast, a.seed, a.dataset_id or "synthetic_demo")
    print(f"SYNTHETIC tutorial dataset written to {out} ({a.n_cases} cases, {a.n_controls} controls)")


if __name__ == "__main__":
    main()
