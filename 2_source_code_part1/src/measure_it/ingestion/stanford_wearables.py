"""Stanford Snyder Lab COVID-19 / Long COVID smartwatch releases (source_id: stanford_longcovid_wearables).

Three public sub-datasets, each namespaced separately (there is no shared
participant key between them):

* ``stanford_covid_mishra2020``  Mishra et al. 2020 Nat Biomed Eng 4:1208 ("phase 1").
  Fitbit minute HR + minute steps (+ sleep) for 118 participants; labels
  (COVID-19 symptom-onset / diagnosis / recovery dates, other illnesses,
  "potential healthy") from Supplementary Data 1. ACUTE COVID-19, no Long COVID label.
* ``stanford_covid_alavi2022``   Alavi et al. 2022 Nat Med 28:175 ("phase 2").
  Fitbit / Apple Watch HR + steps for 2,123 participants; the paper's Source
  Data list 84 COVID-19-positive participants with test date and symptom onset.
  ACUTE COVID-19, no Long COVID label.
* ``stanford_longcovid_uwakwe2025`` Uwakwe et al. 2025 PLOS Digit Health 4:e0001093,
  Stanford Digital Repository druid:cb174pb4851. 126 participants with a
  genuine self-reported Long COVID label (symptoms >= 12 weeks after a
  confirmed diagnosis), demographics, symptoms by phase, and step-masked
  resting HR. The COVID-19 diagnosis date is NOT released, so this subset has
  no days-from-onset alignment.

Run: ``uv run python -m measure_it.ingestion.stanford_wearables``
"""
from __future__ import annotations

import argparse
import json
import re
import warnings
import zipfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from measure_it.config import SEED, TABLES, UNKNOWN, interim_dir, raw_dir, utc_now_iso
from measure_it.download import download_file, load_manifest
from measure_it.provenance import add_provenance
from measure_it.registry import write_registry_entry
from measure_it.store import write_table
from measure_it.wearables import stanford_features as sf

SOURCE_ID = "stanford_longcovid_wearables"
PRODUCER = "measure_it.ingestion.stanford_wearables"
PIPELINE_VERSION = "2026-09-23.3"  # bump to invalidate the per-participant interim cache

GCS = "https://storage.googleapis.com/gbsc-gcp-project-ipop_public"
SDR = "https://stacks.stanford.edu/file/druid:cb174pb4851"

# (url, local filename, subset)
DOWNLOADS = [
    (f"{GCS}/COVID-19/COVID-19-Wearables.zip", "COVID-19-Wearables.zip", "mishra2020"),
    ("https://static-content.springer.com/esm/art%3A10.1038%2Fs41551-020-00640-6/MediaObjects/41551_2020_640_MOESM3_ESM.xlsx",
     "mishra2020_41551_2020_640_MOESM3_ESM.xlsx", "mishra2020"),
    (f"{GCS}/COVID-19-Phase2/COVID-19-Phase2-Wearables.zip", "COVID-19-Phase2-Wearables.zip", "alavi2022"),
    ("https://static-content.springer.com/esm/art%3A10.1038%2Fs41591-021-01593-2/MediaObjects/41591_2021_1593_MOESM3_ESM.xlsx",
     "alavi2022_41591_2021_1593_MOESM3_ESM.xlsx", "alavi2022"),
    (f"{SDR}/Demographic_Data.csv", "longcovid_sdr_Demographic_Data.csv", "uwakwe2025"),
    (f"{SDR}/Symptom_Data.csv", "longcovid_sdr_Symptom_Data.csv", "uwakwe2025"),
    (f"{SDR}/RHR_Data.zip", "longcovid_sdr_RHR_Data.zip", "uwakwe2025"),
    ("https://purl.stanford.edu/cb174pb4851.json", "longcovid_sdr_purl_metadata.json", "uwakwe2025"),
    ("https://journals.plos.org/digitalhealth/article/file?type=supplementary&id=10.1371/journal.pdig.0001093.s009",
     "uwakwe2025_pdig.0001093.s009.xlsx", "uwakwe2025"),
    ("https://journals.plos.org/digitalhealth/article/file?type=supplementary&id=10.1371/journal.pdig.0001093.s010",
     "uwakwe2025_pdig.0001093.s010.xlsx", "uwakwe2025"),
    ("https://journals.plos.org/digitalhealth/article/file?type=supplementary&id=10.1371/journal.pdig.0001093.s011",
     "uwakwe2025_pdig.0001093.s011.xlsx", "uwakwe2025"),
    ("https://storage.googleapis.com/storage/v1/b/gbsc-gcp-project-ipop_public/o?prefix=COVID-19&fields=items(name,size,updated,md5Hash)",
     "gcs_bucket_listing_COVID-19.json", "all"),
]

SUBSETS = {
    "mishra2020": {
        "dataset_id": "stanford_covid_mishra2020",
        "source_name": "Stanford Snyder Lab COVID-19 wearables phase 1 (Mishra et al. 2020, Nat Biomed Eng)",
        "source_version": "GCS object COVID-19/COVID-19-Wearables.zip updated 2020-10-19; labels = Supplementary Data 1 (41551_2020_640_MOESM3_ESM.xlsx)",
        "citation": "Mishra T, et al. Pre-symptomatic detection of COVID-19 from smartwatch data. Nat Biomed Eng. 2020;4:1208-1220. doi:10.1038/s41551-020-00640-6",
        "landing_url": "https://www.nature.com/articles/s41551-020-00640-6",
        "license": "No explicit data licence stated; public download link given in the paper's Data availability statement (Nature Biomedical Engineering 2020).",
    },
    "alavi2022": {
        "dataset_id": "stanford_covid_alavi2022",
        "source_name": "Stanford Snyder Lab COVID-19 wearables phase 2 (Alavi et al. 2022, Nat Med)",
        "source_version": "GCS object COVID-19-Phase2/COVID-19-Phase2-Wearables.zip updated 2021-08-13; labels = Source Data (41591_2021_1593_MOESM3_ESM.xlsx)",
        "citation": "Alavi A, et al. Real-time alerting system for COVID-19 and other stress events using wearable data. Nat Med. 2022;28:175-184. doi:10.1038/s41591-021-01593-2",
        "landing_url": "https://www.nature.com/articles/s41591-021-01593-2",
        "license": "No explicit data licence stated for the zip; article CC BY 4.0; public download link given in the paper's Data availability statement.",
    },
    "uwakwe2025": {
        "dataset_id": "stanford_longcovid_uwakwe2025",
        "source_name": "Snyder Lab Long COVID Study Dataset (Stanford Digital Repository cb174pb4851; Uwakwe et al. 2025, PLOS Digit Health)",
        "source_version": "SDR druid:cb174pb4851 cocina version 3 (deposited 2025-07-12, modified 2026-03-18)",
        "citation": "Uwakwe CK, et al. Longitudinal wearable sensor data enhance precision of Long COVID detection. PLOS Digit Health. 2025;4(11):e0001093. doi:10.1371/journal.pdig.0001093. Data: doi:10.25740/cb174pb4851",
        "landing_url": "https://purl.stanford.edu/cb174pb4851",
        "license": "Open Data Commons Attribution License v1.0 (ODC-By 1.0); SDR use statement: content will not be used to identify individuals.",
    },
}

RHR_DEFINITIONS = {
    "fitbit_minute": "mishra2020_wang: minute-median HR; exclude step minute + next 10 min (WearableDetection rest.min.par=10); no RHR on days with no recorded step",
    "healthkit_sparse": "mishra2020_wang applied to HealthKit samples: minute-median HR; exclude minutes covered by a step interval with steps>0 + next 10 min; no RHR on days with no recorded step",
    "preprocessed_resting_hr": "uwakwe2025_30min_post_step_mask: HR released already step-masked by the authors (30 min after any step)",
}

LONG_COVID_DEFINITION = ("Self-reported symptoms 12 weeks or more after a confirmed (PCR/antigen) COVID-19 diagnosis, "
                         "CDC definition; confirmed by an end-of-study survey (Uwakwe et al. 2025 Methods)")

# Symptom features of the paper's symptoms-only model (RFSM), mapped to Symptom_Data.csv strings.
RFSM_FEATURES = {
    "acute chest pain": ("Chest pain (sharp/aching/burning)", "acute"),
    "acute vomiting": ("Nausea / vomiting", "acute"),
    "post-acute excessive sweating": ("Chills or night sweats", "post_acute"),  # closest available item
    "post-acute memory loss": ("Memory issues", "post_acute"),
    "post-acute brain fog": ("Cognitive dysfunction / slowed thinking / brain fog", "post_acute"),
    "post-acute heart palpitations": ("Palpitations / faster, pounding, or irregular heart beat", "post_acute"),
    "post-acute loss of smell": ("Loss of smell / altered smell", "post_acute"),
}
PHASE_MAP = {
    "Acute COVID-19 [within 4 weeks of COVID-19 detection]": "acute",
    "Post Acute COVID-19 [from 4 to 12 weeks after COVID-19 detection]": "post_acute",
    "Chronic Post COVID-19 Syndrome [12 weeks after COVID-19 detection]": "chronic",
}


# ---------------------------------------------------------------------------
# Download
# ---------------------------------------------------------------------------
def fetch(include_phase2: bool = True) -> dict[str, Path]:
    out = {}
    for url, fn, subset in DOWNLOADS:
        if subset == "alavi2022" and fn.endswith(".zip") and not include_phase2:
            continue
        min_bytes = 100 if fn.endswith((".json", ".csv")) else 1000
        out[fn] = download_file(url, SOURCE_ID, fn, min_bytes=min_bytes, timeout=3600)
    return out


def _retrieved(fn: str) -> str:
    return load_manifest(SOURCE_ID)["files"].get(fn, {}).get("retrieved_at") or utc_now_iso()


# ---------------------------------------------------------------------------
# Labels
# ---------------------------------------------------------------------------
_TS = re.compile(r"Timestamp\('([^']+)'\)")


def parse_timestamp_list(v) -> list[pd.Timestamp]:
    """Supplementary Table 3 stores dates as the text of a Python list of pandas Timestamps."""
    return [pd.Timestamp(x) for x in _TS.findall(str(v))]


def mishra_reference_day(symptoms: list[pd.Timestamp], diagnoses: list[pd.Timestamp]) -> tuple[pd.Timestamp | None, str]:
    """Reference day for a COVID-19 episode in the Mishra release.

    The paper aligns to symptom onset, and to diagnosis when no symptom date
    exists. Supp. Table 3 lists every self-reported symptom date (including
    unrelated earlier illnesses) in no particular order, so the COVID-19 symptom
    onset is taken as the earliest symptom date within 30 days before the
    earliest diagnosis date. Otherwise: the latest symptom date before the
    diagnosis (AIFDJZB, 43 d; matches the paper's Supp. Table 27 day 0), else
    the earliest symptom date; never the list position.
    """
    if diagnoses:
        d0 = min(diagnoses)
        cand = [s for s in symptoms if d0 - pd.Timedelta(days=30) <= s <= d0]
        if cand:
            return min(cand), "symptom_onset"
        before = [s for s in symptoms if s <= d0]
        if before:
            return max(before), "symptom_onset_gt30d_before_diagnosis"
        if symptoms:
            return min(symptoms), "earliest_symptom_date_after_diagnosis"
        return d0, "diagnosis"
    if symptoms:
        return min(symptoms), "earliest_symptom_date"
    return None, "none"


def load_mishra_labels(raw: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    x = raw / "mishra2020_41551_2020_640_MOESM3_ESM.xlsx"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        t3 = pd.read_excel(x, sheet_name="SuppTable3_Fig2a_COVID-19", header=3)
        t29 = pd.read_excel(x, sheet_name="SuppTable29_AlarmFreqComparison", header=2).iloc[:, :3]
    t29.columns = ["label_id", "case_category", "alarm_count_per_30d"]
    t29 = t29.dropna(subset=["label_id"])
    t29["native_id"] = t29.label_id.astype(str).str.replace(r"_\d$", "", regex=True)
    grp = {"COVID positive": "covid19_positive", "Other Illness": "other_illness",
           "potential.healthy": "potential_healthy_no_illness_reported"}
    part = (t29.assign(cohort_group=t29.case_category.map(grp))
            .groupby("native_id").cohort_group.first().reset_index())
    conds = []
    for _, r in t3.dropna(subset=["ParticipantID"]).iterrows():
        lid = str(r.ParticipantID)
        nid = re.sub(r"_Illness_(\d)$", "", lid)
        sym, dia, rec = (parse_timestamp_list(r[c]) for c in ("Symptom_dates", "covid_diagnosis_dates", "recovery_dates"))
        cat = str(r.Category)
        is_covid = cat == "COVID-19"
        # symptom lists are unordered (e.g. AF3J1YC lists 09-03 before 09-02): onset = earliest date, not list position
        ref, basis = mishra_reference_day(sym, dia) if is_covid else ((min(sym), "earliest_symptom_date") if sym else (None, "none"))
        conds.append({
            "native_id": nid,
            "episode_id": lid,
            "condition_type": "infection_episode" if is_covid else "illness_episode",
            "condition_label": "COVID-19 (acute SARS-CoV-2 infection, self-reported diagnosis)" if is_covid
            else f"Non-COVID-19 illness (self-reported){' - ' + cat.split('(')[-1].rstrip(')') if '(' in cat else ''}",
            "label_value": 1,
            "is_long_covid_label": False,
            "symptom_onset_dates": ";".join(str(d.date()) for d in sym) or None,
            "diagnosis_dates": ";".join(str(d.date()) for d in dia) or None,
            "recovery_dates": ";".join(str(d.date()) for d in rec) or None,
            "reference_date": ref.normalize() if ref is not None else pd.NaT,
            "reference_basis": basis,
            "label_source": "Mishra 2020 Supplementary Data 1, sheet SuppTable3_Fig2a_COVID-19",
            "label_definition": "Survey-reported illness with symptom-onset, diagnosis and recovery dates (dates shifted by the authors to obscure PHI, consistently with the wearable files)",
            "symptom_phase": None,
        })
    return part, pd.DataFrame(conds)


def load_alavi_labels(raw: Path) -> pd.DataFrame:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        lab = pd.read_excel(raw / "alavi2022_41591_2021_1593_MOESM3_ESM.xlsx", sheet_name="SourceData_COVID19_Positives")
    lab.columns = ["native_id", "test_date", "symptom_onset", "cohort_timing", "device_reported"]
    lab["test_date"] = pd.to_datetime(lab.test_date.replace("-", np.nan))
    lab["symptom_onset"] = pd.to_datetime(lab.symptom_onset.replace("-", np.nan))
    return lab


def alavi_conditions(lab: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, r in lab.iterrows():
        ref = r.symptom_onset if pd.notna(r.symptom_onset) else r.test_date
        rows.append({
            "native_id": r.native_id, "episode_id": r.native_id, "condition_type": "infection_episode",
            "condition_label": "COVID-19 (acute SARS-CoV-2 infection, positive test)", "label_value": 1,
            "is_long_covid_label": False,
            "symptom_onset_dates": str(r.symptom_onset.date()) if pd.notna(r.symptom_onset) else None,
            "diagnosis_dates": str(r.test_date.date()) if pd.notna(r.test_date) else None,
            "recovery_dates": None,
            "reference_date": pd.Timestamp(ref).normalize() if pd.notna(ref) else pd.NaT,
            "reference_basis": "symptom_onset" if pd.notna(r.symptom_onset) else "test_date_no_symptom_onset_reported",
            "label_source": "Alavi 2022 Source Data, sheet SourceData_COVID19_Positives",
            "label_definition": f"COVID-19 positive test ({r.cohort_timing} cohort); dates shifted by the authors consistently with the wearable files",
            "symptom_phase": None,
        })
    return pd.DataFrame(rows)


def alavi_roster(zip_path: Path) -> pd.DataFrame:
    with zipfile.ZipFile(zip_path) as z:
        names = [n for n in z.namelist() if not n.startswith("__MACOSX") and n.count("/") == 2 and n.endswith(".csv")]
    rows = {}
    for n in names:
        _, pid, fn = n.split("/")
        rows.setdefault(pid, set()).add(fn)
    out = []
    for pid, fns in sorted(rows.items()):
        fit = "Orig_Fitbit_HR.csv" in fns
        out.append({"native_id": pid, "device_class": "fitbit_minute" if fit else "healthkit_sparse",
                    "device": "Fitbit" if fit else "Non-Fitbit (HealthKit export)",
                    "has_hr_file": ("Orig_Fitbit_HR.csv" in fns) or ("Orig_NonFitbit_HR.csv" in fns),
                    "has_steps_file": ("Orig_Fitbit_ST.csv" in fns) or ("Orig_NonFitbit_ST.csv" in fns)})
    return pd.DataFrame(out)


def load_uwakwe(raw: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    demo = pd.read_csv(raw / "longcovid_sdr_Demographic_Data.csv")
    sym = pd.read_csv(raw / "longcovid_sdr_Symptom_Data.csv")
    sym["phase"] = sym.Symptom.str.extract(r"\(choice=(.*)\)\s*$")[0].map(PHASE_MAP)
    sym["symptom"] = sym.Symptom.str.replace(r"\s*\(choice=.*\)\s*$", "", regex=True)
    conds = []
    for _, r in demo.iterrows():
        nid = str(int(r.Participant))
        conds.append({"native_id": nid, "episode_id": nid, "condition_type": "infection_history",
                      "condition_label": "COVID-19 (confirmed prior SARS-CoV-2 infection; date not released)",
                      "label_value": 1, "is_long_covid_label": False, "symptom_onset_dates": None,
                      "diagnosis_dates": None, "recovery_dates": None, "reference_date": pd.NaT,
                      "reference_basis": "not_released",
                      "label_source": "Uwakwe 2025 SDR Demographic_Data.csv (cohort inclusion criterion)",
                      "label_definition": "Confirmed positive COVID-19 test (rapid antigen or PCR) per Methods",
                      "symptom_phase": None})
        conds.append({"native_id": nid, "episode_id": nid, "condition_type": "diagnosis_label",
                      "condition_label": "Long COVID (post-COVID-19 condition)", "label_value": int(r.LongCovid),
                      "is_long_covid_label": True, "symptom_onset_dates": None, "diagnosis_dates": None,
                      "recovery_dates": None, "reference_date": pd.NaT, "reference_basis": "not_released",
                      "label_source": "Uwakwe 2025 SDR Demographic_Data.csv column LongCovid",
                      "label_definition": LONG_COVID_DEFINITION, "symptom_phase": None})
    for _, r in sym.iterrows():
        nid = str(int(r.Participant))
        conds.append({"native_id": nid, "episode_id": nid, "condition_type": "symptom_report",
                      "condition_label": r.symptom, "label_value": int(r.Response), "is_long_covid_label": False,
                      "symptom_onset_dates": None, "diagnosis_dates": None, "recovery_dates": None,
                      "reference_date": pd.NaT, "reference_basis": "not_released",
                      "label_source": "Uwakwe 2025 SDR Symptom_Data.csv",
                      "label_definition": "Symptom reported in the cross-sectional COVID-19 survey for the named phase (acute <4 wk, post-acute 4-12 wk, chronic >=12 wk after detection); only positive responses are released",
                      "symptom_phase": r.phase})
    return demo, sym, pd.DataFrame(conds)


# ---------------------------------------------------------------------------
# Per-participant wearable processing (runs in worker processes)
# ---------------------------------------------------------------------------
_P1_FILE = re.compile(r"^([A-Z0-9]+?)(?:_(\d))?_(hr|steps|sleep)(_longterm)?\.csv$")


def mishra_file_index(zip_path: Path) -> dict[str, dict[str, list[str]]]:
    idx: dict[str, dict[str, list[str]]] = {}
    with zipfile.ZipFile(zip_path) as z:
        for n in z.namelist():
            m = _P1_FILE.match(n.split("/")[-1])
            if m:
                idx.setdefault(m.group(1), {}).setdefault(m.group(3), []).append(n)
    return idx


def wearable_core(hr_min: pd.Series, steps_min: pd.Series | None, device_class: str,
                  run_cusum: bool = True) -> dict:
    """Label-free processing of one participant (cached): daily summary + CuSum daily summary + alarms."""
    th = sf.DEVICE_THRESHOLDS[device_class]
    pre_masked = device_class == "preprocessed_resting_hr"
    daily = sf.daily_summary(hr_min, steps_min, hr_is_resting=pre_masked,
                             min_wear_minutes=th["min_wear_minutes"], min_rhr_minutes=th["min_rhr_minutes"])
    # Without step data the resting-HR definition cannot be applied: report no RHR (hr_mean_all stays)
    rhr_ok = pre_masked or (steps_min is not None and len(steps_min) > 0)
    if len(daily) and not rhr_ok:
        daily[["rhr_mean", "rhr_median", "rhr_night_mean"]] = np.nan
        daily["rhr_minutes"] = 0
        daily["rhr_night_minutes"] = 0
        run_cusum = False
    daily["rhr_computable"] = rhr_ok
    alarms = pd.DataFrame(columns=["alarm_time"])
    evaluable = False
    if len(daily) and run_cusum:
        stats = sf.cusum_online_stats(hr_min, steps_min)
        alarms = sf.cusum_alarms(stats)
        daily = daily.merge(sf.cusum_daily_summary(stats, alarms), on="date", how="left")
        evaluable = stats is not None
    else:
        for c in ("cusum_hours", "cusum_max", "rhr_hourly_res_max", "cusum_alarm"):
            daily[c] = np.nan
    return {"daily": daily, "alarms": alarms, "cusum_evaluable": evaluable, "n_hr_minutes": int(hr_min.size)}


def derive_features(core: dict, reference, device_class: str) -> tuple[pd.DataFrame, dict]:
    """Reference-dependent features (computed every run from the cached label-free core)."""
    daily = core["daily"]
    if daily is None or not len(daily):
        return pd.DataFrame(), {}
    if device_class != "preprocessed_resting_hr":
        # idempotent; also corrects interim caches written before the stepless-day rule existed
        daily = sf.mask_stepless_days(daily)
        # Such older caches get step_stream_present appended last; a fresh core has it right before rhr_computable.
        # Put it there in both cases, so the output schema does not depend on the state of the interim cache
        # (found by the 2026-09-24 cold-interim reproduction: same values, different column position).
        cols = list(daily.columns)
        if "step_stream_present" in cols and "rhr_computable" in cols:
            cols.remove("step_stream_present")
            cols.insert(cols.index("rhr_computable"), "step_stream_present")
            daily = daily[cols]
    daily = sf.add_reference_days(daily, reference)
    base = sf.baseline_stats(daily)
    daily = sf.add_deviations(daily, base)
    feats = {"baseline_days": base.n_days, "baseline_ok": base.ok,
             "baseline_rhr_mean": base.rhr_mean if base.ok else np.nan,
             "baseline_rhr_sd": base.rhr_sd if base.ok else np.nan,
             "baseline_steps_mean": base.steps_mean if base.ok else np.nan}
    feats.update(sf.longitudinal_summary(daily))
    feats.update(sf.recovery_trajectory(daily))
    alarms = core["alarms"]
    ev = bool(core["cusum_evaluable"])
    feats["cusum_evaluable"] = ev
    feats["cusum_alarm_count"] = int(len(alarms)) if ev else np.nan
    feats["cusum_alarms_per_30d"] = float(len(alarms) / max(len(daily), 1) * 30) if ev else np.nan
    for k in ("cusum_alarms_per_30d_pre", "cusum_alarms_per_30d_post", "cusum_infection_alarm_day"):
        feats[k] = np.nan
    if ev and reference is not None and pd.notna(reference):
        ref = pd.Timestamp(reference).normalize()
        t = pd.to_datetime(alarms.alarm_time) if len(alarms) else pd.Series(dtype="datetime64[ns]")
        dd = pd.to_datetime(daily.date)
        pre_days = int((dd < ref - pd.Timedelta(days=14)).sum())
        post_days = int((dd > ref + pd.Timedelta(days=7)).sum())
        feats["cusum_alarms_per_30d_pre"] = float((t < ref - pd.Timedelta(days=14)).sum() / pre_days * 30) if pre_days else np.nan
        feats["cusum_alarms_per_30d_post"] = float((t > ref + pd.Timedelta(days=7)).sum() / post_days * 30) if post_days else np.nan
        feats["cusum_infection_alarm_day"] = sf.select_infection_alarm(alarms, ref)
    return daily, feats


def _cache_path(subset: str, nid: str) -> Path:
    d = interim_dir(SOURCE_ID) / PIPELINE_VERSION / subset
    d.mkdir(parents=True, exist_ok=True)
    return d / f"{nid}.pkl"


def process_mishra(nid: str, files: dict[str, list[str]], zip_path: str) -> dict:
    cp = _cache_path("mishra2020", nid)
    if cp.exists():
        return pd.read_pickle(cp)
    hrs, sts = [], []
    with zipfile.ZipFile(zip_path) as z:
        for n in sorted(files.get("hr", [])):
            h = pd.read_csv(z.open(n), usecols=["datetime", "heartrate"])
            hrs.append(h)
        for n in sorted(files.get("steps", [])):
            s = pd.read_csv(z.open(n), usecols=["datetime", "steps"])
            sts.append(sf.step_minutes_from_rows(s.datetime, s.steps))
    h = pd.concat(hrs, ignore_index=True).drop_duplicates()
    hr_min = sf.minute_median_hr(h.datetime, h.heartrate)
    # short-term and long-term extracts overlap: take the max per minute, never the sum
    steps_min = pd.concat(sts).groupby(level=0).max() if sts else None
    res = {"nid": nid, **wearable_core(hr_min, steps_min, "fitbit_minute")}
    pd.to_pickle(res, cp)
    return res


def _read_member(z: zipfile.ZipFile, name: str, cols: list[str], problems: list[str]) -> pd.DataFrame:
    """Read one CSV member; empty or malformed files become an empty frame and are logged."""
    if name not in z.NameToInfo:
        problems.append(f"missing:{name.split('/')[-1]}")
        return pd.DataFrame(columns=cols)
    try:
        df = pd.read_csv(z.open(name), dtype=str)
    except pd.errors.EmptyDataError:
        problems.append(f"empty_file:{name.split('/')[-1]}")
        return pd.DataFrame(columns=cols)
    miss = [c for c in cols if c not in df.columns]
    if miss:
        head = str(df.columns[0]) if len(df.columns) else ""
        kind = "export_error_text" if head.startswith("BigQuery error") else "missing_columns"
        problems.append(f"{kind}:{name.split('/')[-1]}" + ("" if kind == "export_error_text" else f":{miss}"))
        return pd.DataFrame(columns=cols)
    return df[cols]


def process_alavi(nid: str, device_class: str, zip_path: str) -> dict:
    cp = _cache_path("alavi2022", nid)
    if cp.exists():
        return pd.read_pickle(cp)
    base = f"COVID-19-Phase2-Wearables/{nid}/"
    problems: list[str] = []
    with zipfile.ZipFile(zip_path) as z:
        if device_class == "fitbit_minute":
            h = _read_member(z, base + "Orig_Fitbit_HR.csv", ["datetime", "heartrate"], problems)
            s = _read_member(z, base + "Orig_Fitbit_ST.csv", ["datetime", "steps"], problems)
            steps_min = sf.step_minutes_from_rows(s.datetime, s.steps) if len(s) else None
        else:
            h = _read_member(z, base + "Orig_NonFitbit_HR.csv", ["datetime", "heartrate"], problems)
            s = _read_member(z, base + "Orig_NonFitbit_ST.csv", ["start_datetime", "end_datetime", "steps"], problems)
            steps_min = sf.step_minutes_from_intervals(s.start_datetime, s.end_datetime, s.steps) if len(s) else None
    if steps_min is not None and len(steps_min) == 0:
        steps_min = None
    hr_min = sf.minute_median_hr(h.datetime, h.heartrate)
    if hr_min.empty:
        res = {"nid": nid, "daily": pd.DataFrame(), "alarms": pd.DataFrame(columns=["alarm_time"]),
               "cusum_evaluable": False, "n_hr_minutes": 0, "problems": problems + ["no_valid_hr"],
               "has_steps": steps_min is not None}
        pd.to_pickle(res, cp)
        return res
    # CuSum was developed and validated on Fitbit minute data; not run on sparse HealthKit samples
    res = {"nid": nid, **wearable_core(hr_min, steps_min, device_class, run_cusum=(device_class == "fitbit_minute")),
           "problems": problems, "has_steps": steps_min is not None}
    pd.to_pickle(res, cp)
    return res


def process_uwakwe(nid: str, csv_path: str) -> dict:
    cp = _cache_path("uwakwe2025", nid)
    if cp.exists():
        return pd.read_pickle(cp)
    x = pd.read_csv(csv_path)
    hr_min = sf.minute_median_hr(x.datetime, x.heartrate)
    # RHR already step-masked by the authors; no steps and no infection date -> no CuSum anchoring,
    # but the label-free CuSum alarm rate is still computed on the released resting HR.
    res = {"nid": nid, **wearable_core(hr_min, None, "preprocessed_resting_hr", run_cusum=True)}
    pd.to_pickle(res, cp)
    return res


def _star(args):
    fn, a = args
    return fn(*a)


def _run_parallel(jobs: list, workers: int) -> list[dict]:
    if workers <= 1:
        return [_star(j) for j in jobs]
    with ProcessPoolExecutor(workers) as ex:
        return list(ex.map(_star, jobs, chunksize=1))


# ---------------------------------------------------------------------------
# Reproduction benchmarks
# ---------------------------------------------------------------------------
def reproduce_mishra_cusum(raw: Path, workers: int = 8) -> tuple[list[dict], pd.DataFrame]:
    """CuSum online detection on the short-term phase-1 files vs Supp. Tables 25/27 (Mishra 2020)."""
    x = raw / "mishra2020_41551_2020_640_MOESM3_ESM.xlsx"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        t25 = pd.read_excel(x, sheet_name="SuppTable25_OnlineShortTerm", header=3).iloc[:, :8]
        t27 = pd.read_excel(x, sheet_name="SuppTable27_Fig7f", header=2).iloc[:, :3]
        t3 = pd.read_excel(x, sheet_name="SuppTable3_Fig2a_COVID-19", header=3).set_index("ParticipantID")
    t25.columns = ["pid", "alarm_time", "alarm_hour", "dur_to_max", "total_hours", "max_cusum", "cusum_at_alarm", "pval"]
    t25 = t25.dropna(subset=["pid"])
    t25["alarm_time"] = pd.to_datetime(t25.alarm_time, errors="coerce")
    t27.columns = ["pid", "rhrdiff_day", "cusum_day"]
    t27 = t27.dropna(subset=["pid"])
    zp = str(raw / "COVID-19-Wearables.zip")
    pids = sorted(set(t25.pid.astype(str)))
    jobs = [(_cusum_short_term, (p, zp)) for p in pids]
    res = _run_parallel(jobs, workers)
    ours = pd.concat([r for r in res if len(r)], ignore_index=True) if any(len(r) for r in res) else pd.DataFrame()
    # alarm-level agreement over every participant listed in Supp. Table 25
    n_auth = int(t25.alarm_time.notna().sum())
    key_a = {(p, str(t)) for p, t in zip(t25.pid.astype(str), t25.alarm_time) if pd.notna(t)}
    key_o = {(p, str(t)) for p, t in zip(ours.pid, ours.alarm_time)} if len(ours) else set()
    stat_cols = ["alarm_hour", "dur_to_max", "total_hours", "max_cusum", "cusum_at_alarm"]
    m = t25.dropna(subset=["alarm_time"]).merge(
        ours.rename(columns={"duration_hours_to_max": "dur_to_max"}), on=["pid", "alarm_time"], suffixes=("_paper", "_ours"))
    stats_equal = int(np.all([np.isclose(m[c + "_paper"].astype(float), m[c + "_ours"].astype(float), atol=0.011) for c in stat_cols], axis=0).sum()) if len(m) else 0
    # participant extracts whose alarm lists match exactly (timestamps and run statistics)
    bad_pids = sorted(set(p for p, _ in key_a ^ key_o))
    if len(m):
        ok_stats = np.all([np.isclose(m[c + "_paper"].astype(float), m[c + "_ours"].astype(float), atol=0.011) for c in stat_cols], axis=0)
        bad_pids = sorted(set(bad_pids) | set(m.loc[~ok_stats, "pid"]))
    rows = [{
        "benchmark_id": "mishra2020_cusum_alarm_port",
        "dataset_id": "stanford_covid_mishra2020",
        "paper": SUBSETS["mishra2020"]["citation"],
        "claim": "CuSum online-detection alarms (time, run statistics) listed in Supplementary Table 25 (short-term data, all participant groups)",
        "paper_value": n_auth, "reproduced_value": len(key_a & key_o), "unit": "alarms reproduced exactly (same hour)",
        "n_paper": int(t25.pid.nunique()), "n_reproduced": int(len(pids)),
        "agreement": "agree" if key_a == key_o else "partial",
        "details": f"author alarms {n_auth}; ours {len(key_o)}; identical timestamps {len(key_a & key_o)}; "
                   f"of the matched alarms, {stats_equal} also match alarm hour, hours-to-max, total hours, max CuSum and CuSum-at-alarm (tolerance 0.01); "
                   f"only-in-paper {len(key_a - key_o)}; only-in-ours {len(key_o - key_a)}. "
                   f"Participant extracts reproduced exactly: {len(pids) - len(bad_pids)}/{len(pids)}; discordant: {bad_pids}"
                   + (" (AFHOHOM: the public file covers only 2026-12-31..2027-02-16 (46 HR days) and is absent from COVID-19-Wearables_OLD.zip; the paper's alarm on 2027-01-17 implies the authors analysed a longer extract)" if bad_pids == ["AFHOHOM"] else ""),
        "method": "Python port of WearableDetection R/detection_fn.R (get.stats.fn, stats.fn, cusum.detection.fn; defaults smth.k=10, rest.min=10, base.num=28, resol=60, res.quan=0.9, pval=0.01, max.hour=24, dur.hour=48) run on COVID-19-Wearables.zip <id>_hr.csv/<id>_steps.csv",
        "benchmark_type": "reproduction",
    }]
    # headline: fraction of the 24 COVID-19 cases with >= 28 d pre-onset data that alarm on or before onset
    per = []
    parse = parse_timestamp_list
    for _, r in t27.iterrows():
        p = str(r.pid)
        sym = parse(t3.loc[p, "Symptom_dates"]) if p in t3.index else []
        dia = parse(t3.loc[p, "covid_diagnosis_dates"]) if p in t3.index else []
        ref, basis = mishra_reference_day(sym, dia)
        mine = ours[ours.pid == p] if len(ours) else pd.DataFrame(columns=["alarm_time"])
        theirs = t25[(t25.pid == p) & t25.alarm_time.notna()]
        d_ours = sf.select_infection_alarm(mine, ref)
        d_auth = sf.select_infection_alarm(theirs, ref)
        paper = r.cusum_day
        per.append({"pid": p, "reference_date": ref, "reference_basis": basis,
                    "paper_cusum_day": paper, "ours_cusum_day": d_ours, "rule_on_paper_alarms_day": d_auth})
    per = pd.DataFrame(per)
    paper_num = pd.to_numeric(per.paper_cusum_day, errors="coerce")
    ours_num = per.ours_cusum_day
    paper_le0 = int((paper_num <= 0).sum())
    ours_le0 = int((ours_num <= 0).sum())
    agree_day = int(((paper_num == ours_num) | (paper_num.isna() & ours_num.isna())).sum())
    disagree = per[~((paper_num == ours_num) | (paper_num.isna() & ours_num.isna()))]
    rows.append({
        "benchmark_id": "mishra2020_cusum_presymptomatic_fraction",
        "dataset_id": "stanford_covid_mishra2020",
        "paper": SUBSETS["mishra2020"]["citation"],
        "claim": "62.5% (15/24) of COVID-19 cases with >=28 d of pre-onset data had a CuSum alarm on or before symptom onset (abstract: 63%)",
        "paper_value": round(15 / 24, 4), "reproduced_value": round(ours_le0 / len(per), 4), "unit": "fraction of 24 participants",
        "n_paper": 24, "n_reproduced": int(len(per)),
        "agreement": "agree" if ours_le0 == 15 else "partial",
        "details": f"paper 15/24; recount of paper's Supp. Table 27 gives {paper_le0}/24; ours {ours_le0}/24. "
                   f"Per-participant detection day identical for {agree_day}/24; differing: "
                   + "; ".join(f"{r.pid} paper={r.paper_cusum_day} ours={'miss' if pd.isna(r.ours_cusum_day) else int(r.ours_cusum_day)}" for r in disagree.itertuples()),
        "method": "Infection alarm = alarm whose time falls in the paper's sickness detection window [-14 d, +7 d] around symptom onset (diagnosis when no symptom date), nearest to onset; day = round((alarm - onset)/1 d). The paper does not state its CuSum alarm-to-episode rule; the reference day for multi-episode participants is the earliest symptom date within 30 d before the earliest diagnosis.",
        "benchmark_type": "reproduction",
    })
    # sensitivity to the unstated selection rule (windows x membership x pick)
    vals = []
    for lo in (-14, -21, -28):
        for pick in ("closest", "first"):
            n = 0
            for _, r in per.iterrows():
                mine = ours[ours.pid == r.pid] if len(ours) else pd.DataFrame(columns=["alarm_time"])
                if not len(mine) or pd.isna(r.reference_date):
                    continue
                t = pd.to_datetime(mine.alarm_time)
                ref = pd.Timestamp(r.reference_date)
                w = t[(t >= ref + pd.Timedelta(days=lo)) & (t < ref + pd.Timedelta(days=8))]
                if len(w):
                    off = ((w - ref) / pd.Timedelta(days=1)).to_numpy()
                    d = off[np.argmin(np.abs(off))] if pick == "closest" else off.min()
                    n += int(np.round(d) <= 0)
            vals.append((lo, pick, n))
    rows.append({
        "benchmark_id": "mishra2020_cusum_rule_sensitivity",
        "dataset_id": "stanford_covid_mishra2020",
        "paper": SUBSETS["mishra2020"]["citation"],
        "claim": "Sensitivity of the pre-symptomatic count to the (unstated) alarm-to-episode selection rule",
        "paper_value": 15, "reproduced_value": f"{min(v[2] for v in vals)}-{max(v[2] for v in vals)}",
        "unit": "participants of 24 with alarm day <= 0", "n_paper": 24, "n_reproduced": int(len(per)),
        "agreement": "range_includes_paper" if min(v[2] for v in vals) <= 15 <= max(v[2] for v in vals) else "range_excludes_paper",
        "details": "; ".join(f"window [{lo},+7] {pick}: {n}/24" for lo, pick, n in vals)
                   + ". A matching count is not a matching participant set: the paper counts A3OU183 (alarm at -21.08 d) but calls AYEFCWQ (alarm at -19.04 d) a miss, so no single window rule reproduces the paper's per-participant calls.",
        "method": "Same alarms as mishra2020_cusum_alarm_port; only the episode-matching rule varies",
        "benchmark_type": "sensitivity",
    })
    return rows, per


def _cusum_short_term(pid: str, zip_path: str) -> pd.DataFrame:
    with zipfile.ZipFile(zip_path) as z:
        names = set(z.namelist())
        hn, sn = f"COVID-19-Wearables/{pid}_hr.csv", f"COVID-19-Wearables/{pid}_steps.csv"
        if hn not in names or sn not in names:
            return pd.DataFrame()
        h = pd.read_csv(z.open(hn))
        s = pd.read_csv(z.open(sn))
    al = sf.cusum_alarms(sf.cusum_online_stats(sf.minute_median_hr(h.datetime, h.heartrate),
                                               sf.step_minutes_from_rows(s.datetime, s.steps)))
    al.insert(0, "pid", pid)
    return al


def reproduce_uwakwe(demo: pd.DataFrame, sym: pd.DataFrame, feats: pd.DataFrame, raw: Path,
                     n_splits: int = 1000) -> list[dict]:
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.metrics import average_precision_score, roc_auc_score
    from sklearn.model_selection import StratifiedShuffleSplit

    cit = SUBSETS["uwakwe2025"]["citation"]
    y = demo.set_index("Participant").LongCovid.astype(int)
    rows = []
    # cohort / Table 1
    ct = pd.crosstab(demo.LongCovid, demo.sex)
    vc = pd.crosstab(demo.LongCovid, demo.vaccinated)
    bmi = demo.groupby("LongCovid").bmi.agg(["mean", "std"])
    rows.append({"benchmark_id": "uwakwe2025_cohort_size", "dataset_id": "stanford_longcovid_uwakwe2025", "paper": cit,
                 "claim": "126 participants: 31 with Long COVID, 95 without", "paper_value": "126 (31/95)",
                 "reproduced_value": f"{len(demo)} ({int(y.sum())}/{int((y == 0).sum())})", "unit": "participants",
                 "n_paper": 126, "n_reproduced": int(len(demo)),
                 "agreement": "agree" if (len(demo), int(y.sum())) == (126, 31) else "disagree",
                 "details": "Demographic_Data.csv LongCovid column", "method": "count", "benchmark_type": "reproduction"})
    t1 = (f"female LC {ct.loc[1, 1]}/31 non-LC {ct.loc[0, 1]}/95 (paper 23/31, 42/95); vaccinated-before-infection LC {vc.loc[1, 1]}/31 "
          f"non-LC {vc.loc[0, 1]}/95 (paper 21/31, 81/95); BMI LC {bmi.loc[1, 'mean']:.2f} ({bmi.loc[1, 'std']:.2f}) non-LC "
          f"{bmi.loc[0, 'mean']:.2f} ({bmi.loc[0, 'std']:.2f}) (paper 25.54 (5.95), 27.82 (6.18))")
    ok = (ct.loc[1, 1], ct.loc[0, 1], vc.loc[1, 1], vc.loc[0, 1]) == (23, 42, 21, 81) and \
        abs(bmi.loc[1, "mean"] - 25.54) < 0.01 and abs(bmi.loc[0, "mean"] - 27.82) < 0.01
    rows.append({"benchmark_id": "uwakwe2025_table1", "dataset_id": "stanford_longcovid_uwakwe2025", "paper": cit,
                 "claim": "Table 1 sex, vaccination and BMI by Long COVID status", "paper_value": "see details",
                 "reproduced_value": "see details", "unit": "counts / mean (SD)", "n_paper": 126, "n_reproduced": int(len(demo)),
                 "agreement": "agree" if ok else "disagree", "details": t1 + "; sex coded 1=female (inferred from the exact Table 1 match)",
                 "method": "crosstab of Demographic_Data.csv", "benchmark_type": "reproduction"})
    # symptoms-only random forest (RFSM)
    X = pd.DataFrame(0, index=y.index, columns=list(RFSM_FEATURES))
    for name, (item, phase) in RFSM_FEATURES.items():
        hit = sym[(sym.symptom == item) & (sym.phase == phase)].Participant.unique()
        X.loc[X.index.isin(hit), name] = 1

    def boot(Xm: pd.DataFrame, seed_offset: int) -> tuple[np.ndarray, np.ndarray]:
        sss = StratifiedShuffleSplit(n_splits=n_splits, test_size=0.3, random_state=SEED + seed_offset)
        auc, ap = [], []
        for tr, te in sss.split(Xm, y):
            clf = RandomForestClassifier(n_estimators=200, min_samples_leaf=1, random_state=SEED, n_jobs=4)
            clf.fit(Xm.iloc[tr], y.iloc[tr])
            p = clf.predict_proba(Xm.iloc[te])[:, 1]
            auc.append(roc_auc_score(y.iloc[te], p))
            ap.append(average_precision_score(y.iloc[te], p))
        return np.array(auc), np.array(ap)

    auc, ap = boot(X, 0)

    def fmt(a):
        m, s = a.mean(), a.std(ddof=1)
        return f"{m:.3f} (95% CI of mean {m - 1.96 * s / np.sqrt(a.size):.3f}-{m + 1.96 * s / np.sqrt(a.size):.3f}; 2.5-97.5 pct {np.quantile(a, .025):.3f}-{np.quantile(a, .975):.3f})"
    prev = {c: int(X[c].sum()) for c in X}
    for metric, arr, pv in (("ROC-AUC", auc, 0.907), ("PR-AUC", ap, 0.815)):
        rows.append({"benchmark_id": f"uwakwe2025_rfsm_{metric.lower().replace('-', '_')}", "dataset_id": "stanford_longcovid_uwakwe2025",
                     "paper": cit, "claim": f"Symptoms-only random forest (RFSM) mean test {metric} {pv} over 1,000 bootstrapped 70/30 trials",
                     "paper_value": pv, "reproduced_value": round(float(arr.mean()), 4), "unit": metric,
                     "n_paper": 126, "n_reproduced": int(len(X)),
                     "agreement": "agree" if abs(arr.mean() - pv) <= 0.03 else ("partial" if abs(arr.mean() - pv) <= 0.07 else "disagree"),
                     "details": f"ours {fmt(arr)}; feature prevalence {prev}; 'post-acute excessive sweating' mapped to 'Chills or night sweats' (no sweating item in Symptom_Data.csv)",
                     "method": f"7 RFSM features named in the paper; sklearn RandomForest(200 trees, defaults); {n_splits} stratified 70/30 splits (seed {SEED}). The paper's tuned hyper-parameters and its rule of always keeping its 40% feature-selection subset in training are not reproducible from the release.",
                     "benchmark_type": "reproduction"})
    rows.append({"benchmark_id": "uwakwe2025_rfhrm_rfcm", "dataset_id": "stanford_longcovid_uwakwe2025", "paper": cit,
                 "claim": "HR-only model ROC-AUC 0.748 / PR-AUC 0.462; combined HR+symptom model ROC-AUC 0.951 / PR-AUC 0.859",
                 "paper_value": "0.748 / 0.951", "reproduced_value": UNKNOWN, "unit": "ROC-AUC", "n_paper": 126, "n_reproduced": 0,
                 "agreement": "not_reproducible",
                 "details": "The HR features come from each participant's earliest viable 4-week window after the COVID-19 diagnosis. The diagnosis date is not in the SDR release. The paper's S2 File (windows) and S3 File (202 features) use a participant numbering (1-126, sorted by window start) that does not map onto the SDR ids (0-125): Spearman rho between whole-record mean RR and the published Mean NN = 0.056 (offset +1), 0.008 (offset 0), permutation 97.5% null 0.17; nearest-window feature matching gave non-unique, inexact matches (112/126 unique).",
                 "method": "attempted linkage checks only", "benchmark_type": "reproduction"})
    # exploratory: demographic-only baseline (age bin + sex) under the same resampling
    age_ord = demo.set_index("Participant")["age range"].map({"20-29": 0, "30-39": 1, "40-49": 2, "50-59": 3, "60-69": 4, "70-79": 5})
    Xd = pd.DataFrame({"age_bin": age_ord, "female": demo.set_index("Participant").sex}).reindex(y.index)
    dauc, dap = boot(Xd, 1)
    rows.append({"benchmark_id": "uwakwe2025_exploratory_demographic_rf", "dataset_id": "stanford_longcovid_uwakwe2025", "paper": cit,
                 "claim": "Exploratory (NOT in the paper): demographic-only random forest (10-year age bin + sex) for the Long COVID label",
                 "paper_value": UNKNOWN, "reproduced_value": round(float(dauc.mean()), 4), "unit": "ROC-AUC",
                 "n_paper": 0, "n_reproduced": int(len(Xd)), "agreement": "exploratory",
                 "details": f"ROC-AUC {fmt(dauc)}; PR-AUC {fmt(dap)}; prevalence 31/126 = 0.246 is the PR-AUC chance level",
                 "method": "same resampling as the RFSM row", "benchmark_type": "exploratory_not_reproduction"})
    # exploratory: whole-record resting-HR features (NOT the paper's post-diagnosis window) vs label
    f = feats.set_index("native_id")
    f.index = f.index.astype(int)
    cand = ["rhr_mean_record", "rhr_sd_daily", "rhr_iqr_daily", "rhr_rmssd_daily", "rhr_night_mean_record",
            "rhr_trend_bpm_per_30d", "cusum_alarms_per_30d"]
    rng = np.random.default_rng(SEED)
    for c in cand:
        if c not in f:
            continue
        v = f[c].reindex(y.index)
        ok = v.notna()
        if ok.sum() < 20:
            continue
        a = roc_auc_score(y[ok], v[ok])
        null = np.array([roc_auc_score(rng.permutation(y[ok].values), v[ok]) for _ in range(1000)])
        p = (np.sum(np.abs(null - 0.5) >= abs(a - 0.5)) + 1) / (null.size + 1)
        rows.append({"benchmark_id": f"uwakwe2025_exploratory_auroc_{c}", "dataset_id": "stanford_longcovid_uwakwe2025",
                     "paper": cit, "claim": f"Exploratory (NOT a reproduction): univariate AUROC of whole-record {c} for the Long COVID label",
                     "paper_value": UNKNOWN, "reproduced_value": round(float(a), 4), "unit": "AUROC (0.5 = chance)",
                     "n_paper": 0, "n_reproduced": int(ok.sum()), "agreement": "exploratory",
                     "details": f"two-sided label-permutation p = {p:.4f} (1,000 permutations); record-relative feature, not aligned to infection; no multiple-testing correction across {len(cand)} features",
                     "method": "roc_auc_score on participant-level feature", "benchmark_type": "exploratory_not_reproduction"})
    return rows


# ---------------------------------------------------------------------------
# Assemble tables
# ---------------------------------------------------------------------------
def _prov(df: pd.DataFrame, subset: str, evidence_type: str, retrieved_at: str, notes: str) -> pd.DataFrame:
    s = SUBSETS[subset]
    return add_provenance(df, data_layer="person", source_name=s["source_name"], source_version=s["source_version"],
                          retrieved_at=retrieved_at, evidence_type=evidence_type, source_record_id="native_id",
                          source_geographic_resolution="none", evidence_level="", provenance_notes=notes)


def _pid(subset: str, nid) -> str:
    return f"{SUBSETS[subset]['dataset_id']}:{nid}"


def run(phase2_scope: str = "all", workers: int = 12, n_splits: int = 1000, skip_download: bool = False) -> dict:
    raw = raw_dir(SOURCE_ID)
    if not skip_download:
        fetch()
    started = utc_now_iso()
    rt = {k: _retrieved(fn) for k, fn in (("mishra2020", "COVID-19-Wearables.zip"),
                                          ("alavi2022", "COVID-19-Phase2-Wearables.zip"),
                                          ("uwakwe2025", "longcovid_sdr_RHR_Data.zip"))}

    # ---- labels
    m_part, m_cond = load_mishra_labels(raw)
    a_lab = load_alavi_labels(raw)
    a_cond = alavi_conditions(a_lab)
    a_roster = alavi_roster(raw / "COVID-19-Phase2-Wearables.zip")
    u_demo, u_sym, u_cond = load_uwakwe(raw)

    # ---- phase 1 wearables
    zp1 = raw / "COVID-19-Wearables.zip"
    idx = mishra_file_index(zp1)
    m_ref = (m_cond.sort_values("episode_id").groupby("native_id")
             .apply(lambda g: g.loc[g.condition_type == "infection_episode", "reference_date"].min()
                    if (g.condition_type == "infection_episode").any() else g.reference_date.min(), include_groups=False))
    jobs = [(process_mishra, (nid, idx[nid], str(zp1))) for nid in sorted(idx)]
    m_res = _run_parallel(jobs, workers)

    # ---- phase 2 wearables
    zp2 = raw / "COVID-19-Phase2-Wearables.zip"
    pos = set(a_lab.native_id)
    roster = a_roster if phase2_scope == "all" else a_roster[a_roster.native_id.isin(pos)]
    a_ref = a_cond.set_index("native_id").reference_date
    jobs = [(process_alavi, (r.native_id, r.device_class, str(zp2)))
            for r in roster.itertuples() if r.has_hr_file]
    a_res = _run_parallel(jobs, workers)

    # ---- Long COVID (Uwakwe) wearables
    udir = raw / "longcovid_sdr_RHR_Data"
    if not udir.exists() or len(list(udir.glob("*_hr.csv"))) < len(u_demo):
        with zipfile.ZipFile(raw / "longcovid_sdr_RHR_Data.zip") as z:
            z.extractall(udir)
    jobs = [(process_uwakwe, (str(int(n)), str(udir / f"{int(n)}_hr.csv"))) for n in u_demo.Participant
            if (udir / f"{int(n)}_hr.csv").exists()]
    u_res = _run_parallel(jobs, workers)

    # ---- participants
    parts = []
    mp = m_part.set_index("native_id").cohort_group
    for nid in sorted(set(idx) | set(mp.index)):
        parts.append({"subset": "mishra2020", "native_id": nid, "cohort_group": mp.get(nid, UNKNOWN),
                      "device": "Fitbit", "device_class": "fitbit_minute", "has_wearable_file": nid in idx,
                      "n_hr_files": len(idx.get(nid, {}).get("hr", [])), "n_step_files": len(idx.get(nid, {}).get("steps", [])),
                      "has_sleep_file": bool(idx.get(nid, {}).get("sleep"))})
    for r in a_roster.itertuples():
        parts.append({"subset": "alavi2022", "native_id": r.native_id,
                      "cohort_group": "covid19_positive" if r.native_id in pos else "no_covid19_positive_record_published",
                      "device": r.device, "device_class": r.device_class, "has_wearable_file": bool(r.has_hr_file),
                      "n_hr_files": int(r.has_hr_file), "n_step_files": int(r.has_steps_file), "has_sleep_file": False})
    for r in u_demo.itertuples(index=False):
        d = r._asdict()
        parts.append({"subset": "uwakwe2025", "native_id": str(int(d["Participant"])),
                      "cohort_group": "long_covid" if d["LongCovid"] == 1 else "covid19_recovered_no_long_covid",
                      "device": "mixed (Fitbit / Apple Watch / Garmin / other; per-participant device not released)",
                      "device_class": "preprocessed_resting_hr", "has_wearable_file": (udir / f"{int(d['Participant'])}_hr.csv").exists(),
                      "n_hr_files": 1, "n_step_files": 0, "has_sleep_file": False})
    P = pd.DataFrame(parts)
    processed = {("mishra2020", r["nid"]) for r in m_res if len(r["daily"])} | \
                {("alavi2022", r["nid"]) for r in a_res if len(r["daily"])} | \
                {("uwakwe2025", r["nid"]) for r in u_res if len(r["daily"])}
    P["wearable_processed"] = [(s, n) in processed for s, n in zip(P.subset, P.native_id)]
    probs = {("alavi2022", r["nid"]): ";".join(r.get("problems", [])) for r in a_res}
    P["wearable_read_problems"] = [probs.get((s, n), "") for s, n in zip(P.subset, P.native_id)]
    a_steps = {r["nid"]: bool(r.get("has_steps", False)) for r in a_res}
    P["has_step_data"] = [a_steps.get(n, False) if s == "alavi2022" else (s == "mishra2020" and bool(idx.get(n, {}).get("steps")))
                          for s, n in zip(P.subset, P.native_id)]
    demo = u_demo.rename(columns={"Participant": "native_id"}).copy()
    demo["native_id"] = demo.native_id.astype(int).astype(str)
    demo["sex"] = demo.sex.map({1: "female", 0: "male"})
    demo = demo.rename(columns={"age range": "age_range", "LongCovid": "long_covid_label", "pacific islander": "pacific_islander",
                                "other race": "other_race", "cardiovascular disease": "cardiovascular_disease",
                                "ear nose throat": "ear_nose_throat", "none": "no_comorbidity", "other": "other_condition",
                                "height": "height_in", "weight": "weight_lb"})
    demo.columns = [c.replace(" ", "_") for c in demo.columns]
    P = P.merge(demo.assign(subset="uwakwe2025"), on=["subset", "native_id"], how="left")
    P["sex"] = P["sex"].fillna(UNKNOWN)
    P["age_range"] = P["age_range"].fillna(UNKNOWN)
    P["dataset_id"] = P.subset.map(lambda s: SUBSETS[s]["dataset_id"])
    P["participant_id"] = [_pid(s, n) for s, n in zip(P.subset, P.native_id)]
    P["dates_shifted_by_source"] = P.subset.map({"mishra2020": "yes (Supp. Table notes: dates shifted to obscure PHI)",
                                                 "alavi2022": "yes (shifted consistently with Source Data)",
                                                 "uwakwe2025": "reported ('double date-shifting'); calendar dates not interpretable"})
    P["has_long_covid_label"] = P.subset.eq("uwakwe2025")
    P_out = []
    for s in SUBSETS:
        part = P[P.subset == s].copy()
        P_out.append(_prov(part, s, "person_self_report", rt[s],
                           "Roster of the public release; cohort_group from the paper's supplementary label tables. Person-level; no geography."))
    participants = pd.concat(P_out, ignore_index=True)
    lead = ["participant_id", "dataset_id", "native_id", "cohort_group", "has_long_covid_label", "device", "device_class",
            "has_wearable_file", "wearable_processed"]
    participants = participants[lead + [c for c in participants.columns if c not in lead]]

    # ---- conditions
    C_out = []
    for s, cdf in (("mishra2020", m_cond), ("alavi2022", a_cond), ("uwakwe2025", u_cond)):
        c = cdf.copy()
        c["participant_id"] = [_pid(s, n) for n in c.native_id]
        c["dataset_id"] = SUBSETS[s]["dataset_id"]
        C_out.append(_prov(c, s, "person_self_report", rt[s], "Labels exactly as provided by the source (dates are source-shifted); reference_date is this pipeline's alignment rule"))
    conditions = pd.concat(C_out, ignore_index=True)
    conditions["reference_date"] = pd.to_datetime(conditions.reference_date)
    lead = ["participant_id", "dataset_id", "native_id", "episode_id", "condition_type", "condition_label", "label_value",
            "is_long_covid_label", "symptom_phase", "reference_date", "reference_basis"]
    conditions = conditions[lead + [c for c in conditions.columns if c not in lead]]

    # ---- wearable features + daily
    F_out, D_out = [], []
    for s, res, dclass in (("mishra2020", m_res, None), ("alavi2022", a_res, None), ("uwakwe2025", u_res, "preprocessed_resting_hr")):
        feats, days = [], []
        dc_map = dict(zip(a_roster.native_id, a_roster.device_class)) if s == "alavi2022" else {}
        refs = m_ref if s == "mishra2020" else (a_ref if s == "alavi2022" else pd.Series(dtype="datetime64[ns]"))
        for r in res:
            if not len(r["daily"]):
                continue
            nid = r["nid"]
            dcl = dclass or dc_map.get(nid, "fitbit_minute")
            ref = refs.get(nid, pd.NaT) if len(refs) else pd.NaT
            d, fe = derive_features(r, ref if pd.notna(ref) else None, dcl)
            rdef = RHR_DEFINITIONS[dcl] if bool(r["daily"]["rhr_computable"].any()) else "not computable: no usable step data in the release"
            f = {"native_id": nid, "device_class": dcl, "rhr_definition": rdef,
                 "reference_date": ref, "n_hr_minutes": r.get("n_hr_minutes")}
            f.update(fe)
            feats.append(f)
            d = d.copy()
            d.insert(0, "native_id", nid)
            d["device_class"] = dcl
            days.append(d)
        if not feats:
            continue
        F = pd.DataFrame(feats)
        F["participant_id"] = [_pid(s, n) for n in F.native_id]
        F["dataset_id"] = SUBSETS[s]["dataset_id"]
        if s == "uwakwe2025":
            F["trajectory_label"] = "whole-record resting-HR features; Long COVID label available; NOT aligned to infection (diagnosis date not released)"
        else:
            # a participant without any published illness/infection date has no trajectory, only whole-record features
            F["trajectory_label"] = np.where(
                F.reference_date.notna(),
                "post-infection recovery trajectory (acute COVID-19 / other acute illness); not a Long COVID label",
                "whole-record features only; no published infection/illness date (no trajectory); not a Long COVID label")
        F_out.append(_prov(F, s, "person_derived_feature", rt[s],
                           "Derived by measure_it.wearables.stanford_features; baseline window days -42..-15 before reference; acute peak = max RHR z in days -7..+14; recovery (only if peak z > 2) = first 3 consecutive observed days with RHR z <= 1 from max(peak,0), right-censored at last follow-up (<= day 90)"))
        D = pd.concat(days, ignore_index=True)
        D["participant_id"] = [_pid(s, n) for n in D.native_id]
        D["dataset_id"] = SUBSETS[s]["dataset_id"]
        D_out.append(_prov(D, s, "person_derived_feature", rt[s],
                           "participant x calendar day (source-shifted dates); RHR per the row's device_class definition; deviations only where a pre-reference baseline exists"))
    features = pd.concat(F_out, ignore_index=True)
    daily = pd.concat(D_out, ignore_index=True)
    daily["date"] = pd.to_datetime(daily["date"])
    for c in ("cusum_alarm", "valid_wear_day", "is_baseline_day", "rhr_computable", "step_stream_present"):
        if c in daily:
            daily[c] = daily[c].astype("boolean")
    for c in ("baseline_ok", "cusum_evaluable", "any_elevation_pre_onset_14d", "rhr_acute_elevation"):
        if c in features:
            features[c] = features[c].astype("boolean")
    lead = ["participant_id", "dataset_id", "native_id", "date", "days_from_onset", "device_class"]
    daily = daily[lead + [c for c in daily.columns if c not in lead]]
    features["reference_date"] = pd.to_datetime(features["reference_date"])
    lead = ["participant_id", "dataset_id", "native_id", "device_class", "rhr_definition", "trajectory_label", "reference_date"]
    features = features[lead + [c for c in features.columns if c not in lead]]

    # ---- reproduction benchmarks
    bench, per = reproduce_mishra_cusum(raw, workers=workers)
    ufe = features[features.dataset_id == SUBSETS["uwakwe2025"]["dataset_id"]]
    bench += reproduce_uwakwe(u_demo, u_sym, ufe, raw, n_splits=n_splits)
    bench += [{"benchmark_id": "alavi2022_nightsignal_alerts", "dataset_id": "stanford_covid_alavi2022",
               "paper": SUBSETS["alavi2022"]["citation"], "claim": "Real-time alerts in 67/84 (80%) infected participants (NightSignal)",
               "paper_value": round(67 / 84, 4), "reproduced_value": UNKNOWN, "unit": "fraction", "n_paper": 84, "n_reproduced": 0,
               "agreement": "not_attempted", "details": "NightSignal algorithm (github.com/StanfordBioinformatics/wearable-infection) not ported in this alpha; labels and wearable data are ingested",
               "method": "-", "benchmark_type": "reproduction"}]
    B = pd.DataFrame(bench)
    B["retrieved_at"] = started
    TABLES.mkdir(parents=True, exist_ok=True)
    B.to_csv(TABLES / "stanford_reproduction_benchmark.csv", index=False)
    per.to_csv(raw / "mishra2020_cusum_per_participant_reproduction.csv", index=False)

    # ---- write tables
    write_table(participants, "participants__stanford_covid", producer=PRODUCER,
                description="Roster of the three Stanford Snyder Lab COVID-19/Long COVID wearable releases (dataset-prefixed ids)")
    write_table(conditions, "participant_conditions__stanford_covid", producer=PRODUCER,
                description="Infection/illness episodes with source dates, Long COVID label (Uwakwe 2025 only), symptom reports by phase")
    write_table(features, "participant_wearable_features__stanford_covid", producer=PRODUCER,
                description="One row per participant with processed wearable data: baseline, deviation, post-infection recovery trajectory, longitudinal variability, CuSum alarm rates")
    write_table(daily, "participant_wearable_daily__stanford_covid", producer=PRODUCER,
                description="participant x day: RHR, steps, deviations from pre-reference baseline, days_from_onset, CuSum daily summary")
    summary = {"participants": len(participants), "conditions": len(conditions), "features": len(features),
               "daily": len(daily), "benchmark_rows": len(B), "phase2_scope": phase2_scope}
    write_audit(participants, conditions, features, daily, B, per, summary)
    write_registry(participants, conditions, features, daily, phase2_scope)
    return summary


# ---------------------------------------------------------------------------
# Audit + registry
# ---------------------------------------------------------------------------
def _miss(df: pd.DataFrame, cols: list[str]) -> str:
    lines = ["| variable | n rows | missing | % missing |", "|---|---|---|---|"]
    for c in cols:
        if c in df:
            n = len(df)
            k = int(df[c].isna().sum())
            lines.append(f"| {c} | {n} | {k} | {100 * k / n:.1f} |" if n else f"| {c} | 0 | 0 | - |")
    return "\n".join(lines)


def write_audit(participants, conditions, features, daily, bench, per, summary) -> Path:
    raw = raw_dir(SOURCE_ID)
    man = load_manifest(SOURCE_ID)["files"]
    by = participants.groupby("dataset_id")
    cnt = lambda d, q=None: int(len(participants[(participants.dataset_id == d) & (q if q is not None else True)]))  # noqa: E731
    lines = []
    A = lines.append
    A("# DATA AUDIT — Stanford Snyder Lab COVID-19 / Long COVID wearable releases\n")
    A(f"_Generated by `{PRODUCER}` (pipeline {PIPELINE_VERSION}) at {utc_now_iso()}. Every number below is computed from the retrieved files._\n")
    A("| Field | Value |\n|---|---|")
    A(f"| source_id | {SOURCE_ID} |")
    A("| Source (dataset/API name, exact files/endpoints) | (1) `COVID-19/COVID-19-Wearables.zip` (Mishra 2020) + Nat Biomed Eng Supplementary Data 1; (2) `COVID-19-Phase2/COVID-19-Phase2-Wearables.zip` (Alavi 2022) + Nat Med Source Data; (3) Stanford Digital Repository druid:cb174pb4851 `Demographic_Data.csv`, `Symptom_Data.csv`, `RHR_Data.zip` (Uwakwe 2025) + PLOS Digit Health S2-S4 Files. Bucket: gs://gbsc-gcp-project-ipop_public (public, listable) |")
    A("| Publishing organization | Stanford University School of Medicine, Department of Genetics (Snyder Lab) |")
    A(f"| Retrieval date (UTC) | {min(v['retrieved_at'] for v in man.values())} to {max(v['retrieved_at'] for v in man.values())} |")
    A("| Source version / release | " + "; ".join(f"{v['dataset_id']}: {v['source_version']}" for v in SUBSETS.values()) + " |")
    A("| Source update date / cadence | Static releases: phase 1 zip last-modified 2020-10-19; phase 2 zip 2021-08-13; SDR deposit 2025-07-12, record modified 2026-03-18. No update cadence. |")
    A("| License / access conditions | " + "; ".join(f"{v['dataset_id']}: {v['license']}" for v in SUBSETS.values()) + ". No account, DUA or click-through was needed for any file. |")
    A("| Unit of observation | participant x minute (raw HR/steps); participant x day and participant (processed) |")
    A("| Sample size (actual, as ingested) | " + "; ".join(f"{d}: {len(g)} participants ({int(g.wearable_processed.sum())} with processed wearable data)" for d, g in by) + " |")
    A("| Geography (resolution, vintage) | none. No location fields in any file; participants were recruited across the U.S. and elsewhere via social media. |")
    A("| Person-level? | yes |\n| Geographic? | no |\n| Omics? | no |\n| Wearable? | yes (heart rate, steps; phase-1 sleep files present but not processed) |")
    A("| True participant linkage across modalities? | Within each sub-dataset, wearable files and labels share the same native id. Checked: Mishra 118/118 wearable ids appear in Supp. Table 29; Alavi 84/84 positives appear in the zip; Uwakwe 126/126 RHR files match Demographic_Data ids. **No linkage across the three sub-datasets** (different id schemes, different cohorts/consents). **Not linkable: the Uwakwe paper's S2/S3 files (per-participant windows and features)** use a different participant numbering (see Linkage strategy). |\n")

    A("## Discovery (what exists, and what is a genuine Long COVID label)\n")
    mg = participants[participants.dataset_id == "stanford_covid_mishra2020"].cohort_group.value_counts()
    me = conditions[(conditions.dataset_id == "stanford_covid_mishra2020") & (conditions.condition_type == "illness_episode")]
    ac = conditions[conditions.dataset_id == "stanford_covid_alavi2022"]
    ad = participants[participants.dataset_id == "stanford_covid_alavi2022"].device_class.value_counts()
    ul = conditions[conditions.is_long_covid_label.astype(bool)]
    A("| sub-dataset | URL | participants in release | labels | genuine Long COVID label? |\n|---|---|---|---|---|")
    A(f"| Mishra 2020 phase 1 | {GCS}/COVID-19/COVID-19-Wearables.zip | {cnt('stanford_covid_mishra2020')} | COVID-19 ({mg.get('covid19_positive', 0)}), non-COVID illness ({me.native_id.nunique()} people / {len(me)} episodes), 'potential healthy' ({mg.get('potential_healthy_no_illness_reported', 0)}) with symptom-onset / diagnosis / recovery dates | **No** — acute/post-acute COVID-19 trajectories |")
    A(f"| Alavi 2022 phase 2 | {GCS}/COVID-19-Phase2/COVID-19-Phase2-Wearables.zip | {cnt('stanford_covid_alavi2022')} ({ad.get('fitbit_minute', 0)} Fitbit, {ad.get('healthkit_sparse', 0)} non-Fitbit) | {len(ac)} COVID-19 positives with test date (symptom onset for {int((ac.reference_basis == 'symptom_onset').sum())}); the remaining participants carry no published label (not verified negatives) | **No** — acute COVID-19 |")
    A(f"| Uwakwe 2025 Long COVID (SDR cb174pb4851) | https://purl.stanford.edu/cb174pb4851 | {cnt('stanford_longcovid_uwakwe2025')} | LongCovid 0/1 ({int(ul.label_value.sum())} positive of {len(ul)}), demographics, symptoms by phase | **Yes** — self-reported symptoms >=12 wk after confirmed COVID-19 (CDC definition) |")
    A("\nOther objects in the same public bucket were checked and are not COVID/Long COVID wearable data: `PHD/PHD-paper-cohort-example-data.zip` (insulin-resistance cohort, 5-min HR/steps/sleep with IRvsIS label), `NBME-22-1010C/NBME-22-1010C.zip` (CGM glucose), `HMP/`, `Exercise/`, `Lipidomics/` (iPOP omics), `metabolic_subphenotype_db/`. `COVID-19/COVID-19-Wearables_OLD.zip` is a superseded version of the phase-1 zip and was not ingested.\n")

    A("## Files / endpoints retrieved\n")
    A("| file | url | bytes | sha256 | retrieved_at | last_modified |\n|---|---|---|---|---|---|")
    for fn, v in sorted(man.items()):
        A(f"| {fn} | {v['url']} | {v['bytes']:,} | `{v['sha256'][:16]}…` | {v['retrieved_at']} | {v.get('last_modified') or ''} |")
    A("\nSDR md5 checksums published in the purl metadata matched the downloads (RHR_Data.zip 8282b73b…, Symptom_Data.csv 1d4074f6…, Demographic_Data.csv d5603c76…). The PMC copy of Mishra Supplementary Data 1 returned HTTP 403 to scripted download; the identical Springer ESM file (41551_2020_640_MOESM3_ESM.xlsx) was used.\n")

    A("## Key variables (column headers as read from the files)\n")
    A("* **Mishra phase 1 zip** (`COVID-19-Wearables/<ID>[_n]_{hr,steps,sleep}[_longterm].csv`): `\"user\",\"datetime\",\"heartrate\"` (bpm; minute or 5-15 s samples); `\"user\",\"datetime\",\"steps\"` (per minute, zeros included); `\"user\",\"datetime\",\"stage_duration\",\"stage\"` (sleep; not processed). 280 CSV files (281 zip entries incl. the folder); 118 ids with both HR and steps; 4 ids also have `_longterm` files; AA0HAI1 has 3 per-episode extracts (_1,_2,_3).")
    A("* **Mishra Supplementary Data 1**: `SuppTable3_Fig2a_COVID-19` columns `ParticipantID, Category, Symptom_dates, covid_diagnosis_dates, recovery_dates` (each date cell is the text of a Python list of Timestamps; note 'Dates have been shifted to obscure PHI'); `SuppTable29_AlarmFreqComparison` `ParticipantID, caseCategory, alarmCount` (COVID positive / Other Illness / potential.healthy); `SuppTable25_OnlineShortTerm` and `SuppTable27_Fig7f` used for the benchmark.")
    A("* **Alavi phase 2 zip** (`COVID-19-Phase2-Wearables/<PID>/…`): `Orig_Fitbit_HR.csv` `datetime,heartrate` (5-15 s samples; first line blank); `Orig_Fitbit_ST.csv` `datetime,steps` (only non-zero minutes); `Orig_NonFitbit_HR.csv` `device,datetime,heartrate` (HealthKit; sparse samples); `Orig_NonFitbit_ST.csv` `device,start_datetime,end_datetime,steps` (intervals). participant folders counted in the Discovery table.")
    A("* **Alavi Source Data** sheet `SourceData_COVID19_Positives`: `Participant ID, COVID-19 Test Date, COVID-19 Symptom Onset ('-' = none), Retrospective / Prospective, Device`.")
    A("* **Uwakwe SDR**: `Demographic_Data.csv` `Participant, LongCovid, age range, sex, white, black, native, asian, pacific islander, hispanic, other race, unreported, height, weight, bmi, diabetes, allergy, cardiovascular disease, hypertension, cholesterol, pulmonary, gastrointestinal, hematologic, neurologic, psychiatric, endocrine, musculoskeletal, ear nose throat, genitourinary, dermatologic, lyme, other, none, vaccinated`; `Symptom_Data.csv` `Participant, LongCovid, Symptom, Response` (Symptom = '<item> (choice=<phase>)', Response always 1); `RHR_Data.zip` → `<n>_hr.csv` `datetime,heartrate` (minute averages; already step-masked resting HR, 30-min post-step mask per the paper's Methods — measured: on days with data about 51% of night minutes (00-06 h) vs 11% of midday minutes (11-14 h) are present, i.e. each night hour holds about 8% of all released minutes vs about 1.6% for a midday hour). **No COVID-19 diagnosis/detection date column exists in any SDR file.**")
    A("\nDerived variables (participant_wearable_daily / _features):\n")
    A("| variable | meaning | units |\n|---|---|---|")
    for v, m, u in [
        ("rhr_mean", "mean HR over resting minutes of the day (definition in rhr_definition / module docstring); NaN on days whose step stream recorded no step (step_stream_present = False), where the step filter cannot identify resting minutes", "bpm"),
        ("step_stream_present", "the day has at least one recorded step (phase 1/2 only; NA for the pre-masked Uwakwe HR)", "bool"),
        ("rhr_night_mean", "resting HR 00:00-06:59", "bpm"),
        ("steps", "daily step total on valid wear days (NaN otherwise)", "steps/day"),
        ("valid_wear_day", "HR minutes >= device threshold (Fitbit 600, HealthKit 60) or >= 60 resting minutes (Uwakwe)", "bool"),
        ("days_from_onset", "calendar days from the participant's reference day (symptom onset, else diagnosis/test)", "days"),
        ("rhr_z / rhr_residual", "(rhr_mean - baseline mean)/baseline SD; baseline = days -42..-15, >=7 RHR days", "SD / bpm"),
        ("steps_z / steps_residual", "same for daily steps", "SD / steps"),
        ("cusum_max / cusum_alarm", "daily max of the Mishra/Wang CuSum statistic; alarm flag (Fitbit data only)", "-"),
        ("rhr_acute_peak_z / _day, rhr_acute_elevation", "max daily RHR z in days -7..+14 and whether it exceeds 2", "SD / day / bool"),
        ("days_to_rhr_recovery, rhr_recovery_event", "only if rhr_acute_elevation: first day >= max(peak day, 0) starting 3 consecutive observed days with rhr_z <= 1 (back in the pre-infection band); event=0 = right-censored at last_followup_day (<= day 90)", "days from reference"),
        ("rhr_mean_z_pre14 / _post0_14 / _post15_90, steps_mean_z_*", "window means of daily z relative to the pre-infection baseline", "SD"),
        ("rhr_elevated_days_post", "number of days with z>2 in days 0..90", "days"),
        ("rhr_sd_daily, rhr_rmssd_daily, rhr_iqr_daily, rhr_trend_bpm_per_30d", "whole-record longitudinal variability of daily RHR", "bpm"),
        ("cusum_alarms_per_30d(_pre/_post)", "CuSum alarm rate over the record / before day -14 / after day +7", "alarms per 30 d"),
    ]:
        A(f"| {v} | {m} | {u} |")

    A("\n## Missingness (measured)\n")
    A("### participants__stanford_covid\n")
    A("| dataset | participants | with wearable file | processed | cohort groups |\n|---|---|---|---|---|")
    for d, g in by:
        A(f"| {d} | {len(g)} | {int(g.has_wearable_file.sum())} | {int(g.wearable_processed.sum())} | {g.cohort_group.value_counts().to_dict()} |")
    pr = participants.loc[participants.wearable_read_problems.astype(str) != "", "wearable_read_problems"]
    if len(pr):
        kinds = pr.str.split(";").explode().str.split(":").str[0].value_counts().to_dict()
        A(f"\nPhase-2 file problems (measured while reading): {len(pr)} participants affected; by type {kinds}. "
          "Empty (1-byte) files, files that contain a BigQuery export error message instead of data (export_error_text), and repeated in-file header rows are dropped; "
          "participants without any valid HR are not processed; participants without step data get no RHR and no CuSum (rhr_computable = False).\n")
    nostep = participants[(participants.dataset_id != "stanford_longcovid_uwakwe2025") & participants.wearable_processed & ~participants.has_step_data]
    A(f"\nProcessed participants without usable step data (no RHR/CuSum computed): {len(nostep)} ({nostep.cohort_group.value_counts().to_dict()}).\n")
    u = participants[participants.dataset_id == "stanford_longcovid_uwakwe2025"]
    A("\nUwakwe demographics missingness:\n")
    A(_miss(u, ["height_in", "weight_lb", "bmi", "diabetes", "hypertension", "vaccinated"]))
    A("\n### participant_conditions__stanford_covid\n")
    A("| dataset | condition_type | rows | reference_date missing |\n|---|---|---|---|")
    for (d, t), g in conditions.groupby(["dataset_id", "condition_type"]):
        A(f"| {d} | {t} | {len(g)} | {int(g.reference_date.isna().sum())} |")
    rb = conditions[conditions.dataset_id != "stanford_longcovid_uwakwe2025"].groupby(["dataset_id", "condition_type"]).reference_basis.value_counts()
    A("\nreference_basis (alignment rule actually used): " + "; ".join(f"{d}/{t}/{b}: {n}" for (d, t, b), n in rb.items()) + "\n")
    A("\n### participant_wearable_daily__stanford_covid\n")
    for d, g in daily.groupby("dataset_id"):
        nostep = g.step_stream_present.eq(False) & g.rhr_computable.fillna(False).astype(bool) if "step_stream_present" in g else pd.Series(False, index=g.index)
        A(f"\n**{d}** — {len(g):,} participant-days, {g.participant_id.nunique()} participants; valid wear days {int(g.valid_wear_day.fillna(False).sum()):,}; days with days_from_onset {int(g.days_from_onset.notna().sum()):,}; "
          f"days with step data but no recorded step (RHR and steps blanked) {int(nostep.sum()):,} ({int((nostep & g.valid_wear_day.fillna(False).astype(bool)).sum()):,} of them valid wear days, {g.loc[nostep, 'participant_id'].nunique()} participants)\n")
        A(_miss(g, ["rhr_mean", "rhr_night_mean", "steps", "rhr_z", "steps_z", "cusum_max"]))
    A("\n### participant_wearable_features__stanford_covid\n")
    for d, g in features.groupby("dataset_id"):
        A(f"\n**{d}** — {len(g)} participants; reference date {int(g.reference_date.notna().sum())}; baseline_ok {int(g.baseline_ok.fillna(False).sum())}; "
          f"acute RHR elevation (peak z>2) {int(g.rhr_acute_elevation.fillna(False).sum())}; recovery observed {int((g.rhr_recovery_event == 1).sum())}, right-censored {int((g.rhr_recovery_event == 0).sum())}\n")
        A(_miss(g, ["reference_date", "baseline_rhr_mean", "rhr_acute_peak_z", "days_to_rhr_recovery", "rhr_sd_daily", "steps_mean_daily", "cusum_alarms_per_30d"]))

    A("\n## Linkage strategy\n")
    A("* `participant_id = <dataset_id>:<native id>` with dataset_id in {stanford_covid_mishra2020, stanford_covid_alavi2022, stanford_longcovid_uwakwe2025}. Join labels and wearables only within a dataset_id.")
    A("* Mishra: wearable filename id == Supp. Table 3/29 ParticipantID (AA0HAI1_1..3 in the tables = episodes of file id AA0HAI1).")
    A("* Alavi: folder name == Source Data `Participant ID` (84/84 found).")
    A("* Uwakwe: `<n>_hr.csv` == Demographic_Data/Symptom_Data `Participant` (0-125). The paper's S2 File (4-week windows, weeks after diagnosis) and S3 File (202 HR features) use `Participant` 1-126 sorted by window start; tested and REJECTED as a mapping (Spearman rho whole-record mean RR vs published Mean NN = 0.056 at offset +1, 0.008 at offset 0, 0.053 at offset -1; label-permutation null 97.5% = 0.17). A 4-feature nearest-window search found no exact matches (112 distinct targets for 126 participants). The file start is not the diagnosis date either (58 participants would start at week 0 vs 10 in S2; KS p = 1e-14). Hence no infection alignment is possible for this subset.")
    A("* No person in these data can be linked to any geographic or facility row; there are no location fields.\n")

    A("## Limitations and caveats\n")
    A("* **Only the Uwakwe 2025 subset carries a Long COVID label.** Mishra and Alavi are acute COVID-19 cohorts; their features are named 'post-infection recovery trajectory', not Long COVID.")
    A("* Uwakwe: the label is self-reported symptom duration (recall bias); age is released only as 10-year bins; the comorbidity block is missing for 84/126; per-participant device is not released; HR is already masked by the authors' 30-min post-step rule; steps are not released; no infection date, so trajectories cannot be aligned to infection and the HR-only/combined models cannot be reproduced.")
    A("* All calendar dates are shifted by the source (random per-person offsets, some into 2027-2030). Day-of-week is not guaranteed and seasonality is not interpretable; only within-person intervals are meaningful.")
    A("* Selection: highly engaged volunteers with smartwatches (mostly White per Uwakwe Table 1); Mishra analysed Fitbit users only; participants who stop wearing the device when ill produce gaps exactly around onset.")
    A("* Alavi non-positive participants have no published label; they are **not** verified COVID-negative controls.")
    A("* Apple Watch (HealthKit) HR is sparse (about 50-360 minutes per day) so its RHR is noisier and uses lower thresholds; CuSum is only run on Fitbit minute data (the validated input) and on the Uwakwe resting HR (label-free alarm rate only).")
    A("* Baseline (days -42..-15) needs >=7 RHR days; participants whose data start after day -15 have no deviations or recovery values. Recovery is defined only for participants whose daily RHR z exceeded 2 in days -7..+14, and is right-censored when follow-up ends before 3 consecutive in-band days. A single-day z>2 in a 22-day window also occurs by chance, so rhr_acute_elevation is descriptive, not a detection claim.")
    A("* The CuSum pre-symptomatic count depends on an alarm-to-episode matching rule the paper does not state (see benchmark sensitivity row).")
    A("* Step-stream dropouts: on worn days with no recorded step, all-minute HR equals the participant's usual all-minute HR (median difference about -0.1 bpm) and sits 6-7 bpm above their usual RHR, so these are dropouts, not rest. Daily RHR/steps are blanked on those days (counts in Missingness). The CuSum columns follow the authors' published algorithm unchanged, which treats such days as all-resting, so cusum_* can still react to a dropout.")
    A("* Daily RHR does not require a valid wear day: a partially worn day (e.g. night only) with enough resting minutes still yields an RHR, which is biased toward night-time values.\n")

    A("## Reproduction benchmark (results/tables/stanford_reproduction_benchmark.csv)\n")
    A("| benchmark | paper | ours | agreement | details |\n|---|---|---|---|---|")
    for r in bench.itertuples():
        A(f"| {r.benchmark_id} | {r.paper_value} | {r.reproduced_value} | {r.agreement} | {str(r.details)[:400]} |")
    A("\nAdditional development check (not rerun by `run()`): the CuSum port run on the WearableDetection repository's own example input (data/AHYIJDV_hr.csv, data/AHYIJDV_steps.csv) reproduced its published result/CuSum_online_detection.csv exactly (3/3 alarms: 2020-09-25 10:00, 2020-10-13 07:00, 2020-10-22 01:00 with identical run statistics).\n")

    A("## Processed outputs\n")
    for name, df in (("participants__stanford_covid", participants), ("participant_conditions__stanford_covid", conditions),
                     ("participant_wearable_features__stanford_covid", features), ("participant_wearable_daily__stanford_covid", daily)):
        A(f"* `data/processed/{name}.parquet` — {len(df):,} rows")
    A(f"* `results/tables/stanford_reproduction_benchmark.csv` — {len(bench)} rows")
    A("* `data/raw/stanford_longcovid_wearables/mishra2020_cusum_per_participant_reproduction.csv` — per-participant CuSum detection day, paper vs ours\n")
    A("## Reproduce\n")
    A("`uv run python -m measure_it.ingestion.stanford_wearables` (downloads ~6.1 GB on first run; per-participant results cached under data/interim/stanford_longcovid_wearables/). Options: `--phase2-scope positives` to process only the 84 labelled phase-2 participants; `--workers N`.")
    p = raw / "DATA_AUDIT.md"
    p.write_text("\n".join(lines) + "\n")
    return p


def write_registry(participants, conditions, features, daily, phase2_scope) -> Path:
    man = load_manifest(SOURCE_ID)["files"]
    n = participants.groupby("dataset_id").size().to_dict()
    npro = participants.groupby("dataset_id").wearable_processed.sum().astype(int).to_dict()
    entry = {
        "source_id": SOURCE_ID,
        "name": "Stanford Snyder Lab COVID-19 and Long COVID smartwatch datasets (Mishra 2020, Alavi 2022, Uwakwe 2025)",
        "publisher": "Stanford University (Snyder Lab, Department of Genetics)",
        "landing_url": "https://snyderlabs.stanford.edu/long-covid/",
        "access_urls": sorted({v["url"] for v in man.values()}),
        "license": " | ".join(f"{v['dataset_id']}: {v['license']}" for v in SUBSETS.values()),
        "access_conditions": "Open download, no registration or DUA (verified 2026-09-23). SDR terms: do not use to identify individuals.",
        "retrieved_at": max(v["retrieved_at"] for v in man.values()),
        "source_version": " | ".join(f"{v['dataset_id']}: {v['source_version']}" for v in SUBSETS.values()),
        "update_date": "phase1 zip 2020-10-19; phase2 zip 2021-08-13; SDR modified 2026-03-18",
        "data_layer": "person",
        "unit_of_observation": "participant x minute (raw); participant x day; participant",
        "sample_size": {"participants": n, "participants_with_processed_wearables": npro,
                        "participant_days": int(len(daily)),
                        "long_covid_labelled": int(conditions.is_long_covid_label.astype(bool).sum()),
                        "long_covid_positive": int(conditions.loc[conditions.is_long_covid_label.astype(bool), "label_value"].sum()),
                        "acute_covid19_positive_with_dates": {k: int(v) for k, v in conditions[conditions.condition_type == "infection_episode"]
                                                              .groupby("dataset_id").native_id.nunique().items()}},
        "geographic_resolution": "none",
        "person_level": True, "geographic": False, "omics": False, "wearable": True,
        "participant_linkage": "Within each sub-dataset wearable files and labels share native ids (verified); no linkage across sub-datasets; Uwakwe paper supplementary windows/features are on an unlinkable numbering",
        "true_participant_linkage_across_modalities": False,
        "status": "ingested",
        "processed_outputs": ["participants__stanford_covid", "participant_conditions__stanford_covid",
                              "participant_wearable_features__stanford_covid", "participant_wearable_daily__stanford_covid",
                              "results/tables/stanford_reproduction_benchmark.csv"],
        "audit": f"data/raw/{SOURCE_ID}/DATA_AUDIT.md",
        "ingestion_module": PRODUCER,
        "limitations": [
            "Only stanford_longcovid_uwakwe2025 (n=126, 31 positive) has a genuine Long COVID label; Mishra/Alavi are acute COVID-19 cohorts",
            "Uwakwe release has no COVID-19 diagnosis date and no steps: no infection-aligned trajectories; HR-only and combined models not reproducible",
            "All dates are source-shifted; calendar effects not interpretable",
            "Alavi participants without a published positive record are not verified negatives",
            f"phase-2 processing scope for this run: {phase2_scope}",
        ],
        "sub_datasets": {v["dataset_id"]: {"citation": v["citation"], "landing_url": v["landing_url"]} for v in SUBSETS.values()},
        "notes": "Long COVID label = self-reported symptoms >=12 weeks after confirmed COVID-19 (CDC). RHR = Mishra/Wang CuSum-code definition (step minute + 10 min excluded; no RHR on days with no recorded step) except Uwakwe (authors' 30-min mask).",
    }
    return write_registry_entry(entry)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--phase2-scope", choices=["all", "positives"], default="all")
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--n-splits", type=int, default=1000)
    ap.add_argument("--skip-download", action="store_true")
    a = ap.parse_args()
    print(json.dumps(run(phase2_scope=a.phase2_scope, workers=a.workers, n_splits=a.n_splits,
                         skip_download=a.skip_download), indent=2, default=str))


if __name__ == "__main__":
    main()
