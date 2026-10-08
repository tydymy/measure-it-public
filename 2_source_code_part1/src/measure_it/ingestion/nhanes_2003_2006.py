"""NHANES 2003-2004 (C) and 2005-2006 (D): person-level clinical ingestion for a hip-accelerometry replication cohort.

These cycles carry hip-worn ActiGraph AM-7164 accelerometry (PAXRAW_C/D: minute "intensity" counts, plus steps in
2005-2006) — a different device, body site and unit from the 2011-2014 wrist MIMS data. The accelerometry is
processed by ``measure_it.wearables.nhanes0306_features``; this module ingests demographics, the questionnaire,
examination and laboratory components that match the 2011-2014 analysis (CRP exists in these cycles), prescription
medications and the public-use linked mortality files.

Sources (verified 2026-09-24):
  * https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/{2003|2005}/DataFiles/<FILE>.xpt (+ <FILE>.htm documentation)
  * https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/{2003|2005}/DataFiles/PAXRAW_{C|D}.zip (listed on the examination
    data pages; the .xpt path returns HTTP 404)
  * https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality/NHANES_{2003_2004|2005_2006}_MORT_2019_PUBLIC.dat

participant_id = "nhanes0306:<SEQN>" — a namespace distinct from the 2011-2014 "nhanes:" ids (SEQNs do not
overlap across cycles, but the datasets are kept separate). Joins are made ONLY on SEQN.

Outputs: participants__nhanes0306, participant_clinical_features__nhanes0306, participant_labs__nhanes0306,
participant_medications__nhanes0306, participant_mortality__nhanes0306; DATA_AUDIT.md and registry_entry.yaml under
data/raw/nhanes_2003_2006/.

Run: uv run python -m measure_it.ingestion.nhanes_2003_2006
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

from ..config import UNKNOWN, raw_dir, utc_now_iso
from ..download import download_file, load_manifest
from ..provenance import add_provenance
from ..registry import write_registry_entry
from ..store import write_table
from . import nhanes as N

SOURCE_ID = "nhanes_2003_2006"
SOURCE_NAME = "NHANES 2003-2004 (C) and 2005-2006 (D), CDC NCHS"
PRODUCER = "measure_it.ingestion.nhanes_2003_2006"
DATASET_PREFIX = "nhanes0306"
CYCLES = {"C": {"year": 2003, "label": "2003-2004"}, "D": {"year": 2005, "label": "2005-2006"}}
PUBLIC_BASE = N.PUBLIC_BASE
MORT_FILES = {"C": "NHANES_2003_2004_MORT_2019_PUBLIC.dat", "D": "NHANES_2005_2006_MORT_2019_PUBLIC.dat"}
LISTING = "https://wwwn.cdc.gov/nchs/nhanes/search/datapage.aspx?Component={comp}&CycleBeginYear={year}"

# role -> file per cycle (None = not published in that cycle; recorded in the audit)
FILES = {
    "DEMO": {"C": "DEMO_C", "D": "DEMO_D"},
    "BPX": {"C": "BPX_C", "D": "BPX_D"}, "BMX": {"C": "BMX_C", "D": "BMX_D"},
    "MCQ": {"C": "MCQ_C", "D": "MCQ_D"}, "BPQ": {"C": "BPQ_C", "D": "BPQ_D"}, "CDQ": {"C": "CDQ_C", "D": "CDQ_D"},
    "DIQ": {"C": "DIQ_C", "D": "DIQ_D"}, "PFQ": {"C": "PFQ_C", "D": "PFQ_D"}, "HSQ": {"C": "HSQ_C", "D": "HSQ_D"},
    "HUQ": {"C": "HUQ_C", "D": "HUQ_D"}, "SMQ": {"C": "SMQ_C", "D": "SMQ_D"}, "RXQ_RX": {"C": "RXQ_RX_C", "D": "RXQ_RX_D"},
    "SLQ": {"C": None, "D": "SLQ_D"}, "DPQ": {"C": None, "D": "DPQ_D"},
    "CBC": {"C": "L25_C", "D": "CBC_D"}, "BIOPRO": {"C": "L40_C", "D": "BIOPRO_D"}, "GHB": {"C": "L10_C", "D": "GHB_D"},
    "TCHOL": {"C": "L13_C", "D": "TCHOL_D"}, "HDL": {"C": "L13_C", "D": "HDL_D"},
    "VID": {"C": "VID_C", "D": "VID_D"}, "B12": {"C": "L06NB_C", "D": "B12_D"},
    "CRP": {"C": "L11_C", "D": "CRP_D"}, "FERRITIN": {"C": "L06TFR_C", "D": "FERTIN_D"},
    "GLU": {"C": "L10AM_C", "D": "GLU_D"}, "TRIGLY": {"C": "L13AM_C", "D": "TRIGLY_D"},
}
NOT_PUBLISHED = {
    "DPQ_C": "PHQ-9 (DPQ) starts in 2005-2006; 2003-2004 used the CIDI depression module (CIQDEP_C, ages 20-39 only) "
             "-> fatigue, PHQ-9 and the ME/CFS-like proxy are 2005-2006 only",
    "SLQ_C": "sleep questionnaire (SLQ) starts in 2005-2006; not on the 2003 questionnaire listing",
    "INQ_C/D": "income questionnaire INQ starts 2007; the DEMO income-poverty ratio INDFMPIR is used",
    "MCQ082 / MCQ160O": "celiac disease (MCQ082) and COPD (MCQ160O) are not asked in 2003-2006",
    "DLQ": "disability questionnaire (DLQ) starts in 2013-2014",
    "LBXSCK": "creatine kinase is not part of the 2003-2006 standard biochemistry profile (L40_C / BIOPRO_D); the "
              "2011-2014 clinical feature lab_ck_iu_l is therefore all-missing here and dropped by the imputer",
    "VID_C/VID_D": "vitamin D is released as LBDVIDMS (NCHS 2015 LC-MS/MS-equivalent re-release; variable name "
                   "differs from 2011-2014 LBXVIDMS)",
}

# analysis lab variables: harmonised name -> NHANES variable (same in both cycles unless noted)
LAB_WIDE = {
    "LBXWBCSI": "wbc_1000_per_ul", "LBXLYPCT": "lymphocyte_pct", "LBXNEPCT": "neutrophil_pct",
    "LBXEOPCT": "eosinophil_pct", "LBXHGB": "hemoglobin_g_dl", "LBXHCT": "hematocrit_pct", "LBXMCVSI": "mcv_fl",
    "LBXRDW": "rdw_pct", "LBXPLTSI": "platelets_1000_per_ul", "LBXSAL": "albumin_g_dl", "LBXSATSI": "alt_u_l",
    "LBXSASSI": "ast_u_l", "LBXSAPSI": "alp_u_l", "LBXSBU": "bun_mg_dl", "LBXSCR": "creatinine_mg_dl",
    "LBXSCK": "ck_iu_l", "LBXSLDSI": "ldh_u_l", "LBXSUA": "uric_acid_mg_dl", "LBXSGL": "glucose_serum_mg_dl",
    "LBXSIR": "iron_ug_dl", "LBXSNASI": "sodium_mmol_l", "LBXSKSI": "potassium_mmol_l",
    "LBXSC3SI": "bicarbonate_mmol_l", "LBXSTP": "total_protein_g_dl", "LBXSGB": "globulin_g_dl",
    "LBXSCA": "calcium_mg_dl", "LBXSGTSI": "ggt_u_l", "LBXGH": "hba1c_pct", "LBXTC": "total_cholesterol_mg_dl",
    "HDL": "hdl_mg_dl", "LBDVIDMS": "vitamin_d_25oh_nmol_l", "LBDFER": "ferritin_ng_ml", "LBXB12": "vitamin_b12_pg_ml",
    "LBXCRP": "crp_mg_dl", "LBXFER": "ferritin_ng_ml", "LBXGLU": "glucose_fasting_mg_dl",
    "LBXTR": "triglycerides_fasting_mg_dl", "LBDLDL": "ldl_fasting_mg_dl",
}
HDL_VARS = {"C": "LBXHDD", "D": "LBDHDD"}
LAB_ROLE_VARS = {"CBC": ["LBXWBCSI", "LBXLYPCT", "LBXNEPCT", "LBXEOPCT", "LBXHGB", "LBXHCT", "LBXMCVSI", "LBXRDW",
                         "LBXPLTSI"],
                 "BIOPRO": ["LBXSAL", "LBXSATSI", "LBXSASSI", "LBXSAPSI", "LBXSBU", "LBXSCR", "LBXSCK", "LBXSLDSI",
                            "LBXSUA", "LBXSGL", "LBXSIR", "LBXSNASI", "LBXSKSI", "LBXSC3SI", "LBXSTP", "LBXSGB",
                            "LBXSCA", "LBXSGTSI"],
                 "GHB": ["LBXGH"], "TCHOL": ["LBXTC"], "HDL": ["HDL"], "VID": ["LBDVIDMS"], "B12": ["LBXB12"],
                 "CRP": ["LBXCRP"], "FERRITIN": ["LBXFER", "LBDFER"], "GLU": ["LBXGLU"], "TRIGLY": ["LBXTR", "LBDLDL"]}
LAB_SUBSAMPLE = {"FERRITIN": "ferritin measured in a subgroup (see codebook target; 2003-2004 file also holds "
                             "transferrin receptor)",
                 "GLU": "morning fasting subsample (weight WTSAF2YR)", "TRIGLY": "morning fasting subsample (weight WTSAF2YR)"}

# items as in measure_it.ingestion.nhanes (same output names); items absent in a cycle stay NaN
YN_ITEMS = [(c, v, n) for c, v, n in N.YN_ITEMS if c not in ("DLQ", "HIQ", "PAQ")]
RANGE_ITEMS = [(c, v, n, lo, hi) for c, v, n, lo, hi in N.RANGE_ITEMS if c not in ("INQ", "PAQ")]
RACE1 = {1: "Mexican American", 2: "Other Hispanic", 3: "Non-Hispanic White", 4: "Non-Hispanic Black",
         5: "Other race including multi-racial"}


def pid(seqn) -> pd.Series:
    return DATASET_PREFIX + ":" + pd.Series(seqn).astype("int64").astype(str)


def _url(file: str, cyc: str, ext: str = "xpt") -> str:
    return PUBLIC_BASE.format(year=CYCLES[cyc]["year"], file=f"{file}.{ext}")


def fetch_all(log=print) -> None:
    (raw_dir(SOURCE_ID) / "docs" / "listings").mkdir(parents=True, exist_ok=True)
    for cyc, meta in CYCLES.items():
        for comp in ("Laboratory", "Examination", "Questionnaire", "Demographics"):
            download_file(LISTING.format(comp=comp, year=meta["year"]), SOURCE_ID,
                          f"docs/listings/{comp.lower()}_{meta['year']}.html", allow_html=True, min_bytes=20000)
    seen = set()
    for role, per in FILES.items():
        for cyc, f in per.items():
            if not f or f in seen:
                continue
            seen.add(f)
            download_file(_url(f, cyc), SOURCE_ID, f"{f}.xpt", min_bytes=500, max_retries=4)
            download_file(_url(f, cyc, "htm"), SOURCE_ID, f"docs/{f}.htm", allow_html=True, min_bytes=3000,
                          max_retries=4)
    for cyc in CYCLES:
        download_file(_url(f"PAXRAW_{cyc}", cyc, "htm"), SOURCE_ID, f"docs/PAXRAW_{cyc}.htm", allow_html=True,
                      min_bytes=3000)
        download_file(_url(f"PAXRAW_{cyc}", cyc, "zip"), SOURCE_ID, f"PAXRAW_{cyc}.zip", min_bytes=10_000_000,
                      timeout=1200, max_retries=4)
        download_file(N.MORT_BASE.format(file=MORT_FILES[cyc]), SOURCE_ID, MORT_FILES[cyc], min_bytes=1000)
    log("[fetch] all 2003-2006 files present")


def xpt(f: str) -> Path:
    return raw_dir(SOURCE_ID) / f"{f}.xpt"


def codebook(f: str) -> dict:
    comp, _, cyc = f.rpartition("_")
    return N.parse_codebook(comp, cyc, source_id=SOURCE_ID)


def load(role: str, columns: list[str]) -> pd.DataFrame:
    frames = []
    for cyc, f in FILES[role].items():
        if not f:
            continue
        lay = N.xport_layout(xpt(f))
        have = [c for c in columns if c.upper() in lay["by_name"] and c.upper() != "SEQN"]
        df = N.read_xport(xpt(f), ["SEQN"] + have)
        for c in columns:
            if c not in df.columns:
                df[c] = np.nan
        df["cycle_code"] = cyc
        frames.append(df)
    out = pd.concat(frames, ignore_index=True)
    out["SEQN"] = out["SEQN"].astype("int64")
    return out


def cycle_meta(cyc: str) -> dict:
    m = load_manifest(SOURCE_ID)["files"].get(f"DEMO_{cyc}.xpt", {})
    return {"source_version": f"NHANES {CYCLES[cyc]['label']} public release (cycle {cyc}); DEMO_{cyc}.xpt "
                              f"last-modified {m.get('last_modified', UNKNOWN)}",
            "retrieved_at": m.get("retrieved_at", UNKNOWN)}


def write_person(df: pd.DataFrame, name: str, *, evidence_type: str, notes: str, description: str) -> Path:
    parts = []
    for cyc, sub in df.groupby("cycle_code", sort=True):
        meta = cycle_meta(cyc)
        parts.append(add_provenance(sub, data_layer="person", source_name=SOURCE_NAME,
                                    source_version=meta["source_version"], retrieved_at=meta["retrieved_at"],
                                    evidence_type=evidence_type, source_record_id="record_key",
                                    source_geographic_resolution="none", provenance_notes=notes))
    out = pd.concat(parts, ignore_index=True).drop(columns=["record_key"]).rename(columns={"SEQN": "seqn"})
    return write_table(out, name, producer=PRODUCER, description=description)


# --------------------------------------------------------------------------------------------- builders


def build_participants() -> pd.DataFrame:
    d = load("DEMO", ["SDDSRVYR", "RIDSTATR", "RIAGENDR", "RIDAGEYR", "RIDRETH1", "DMDEDUC2", "INDFMPIR", "RIDEXMON",
                      "RIDEXPRG", "SDMVPSU", "SDMVSTRA", "WTINT2YR", "WTMEC2YR"])
    out = pd.DataFrame({
        "participant_id": pid(d["SEQN"]).values, "SEQN": d["SEQN"], "cycle_code": d["cycle_code"],
        "cycle": d["cycle_code"].map(lambda c: CYCLES[c]["label"]), "sddsrvyr": d["SDDSRVYR"].astype("Int64"),
        "mec_examined": d["RIDSTATR"].eq(2.0), "sex": d["RIAGENDR"].map({1.0: "male", 2.0: "female"}),
        "age_years": d["RIDAGEYR"], "age_topcoded_85": d["RIDAGEYR"].eq(85.0),
        "race_ethnicity": d["RIDRETH1"].map(RACE1), "education_adult20plus": d["DMDEDUC2"].map(N.EDUC2),
        "income_poverty_ratio": d["INDFMPIR"], "exam_period": d["RIDEXMON"].map({1.0: "Nov-Apr", 2.0: "May-Oct"}),
        "pregnant_at_exam": d["RIDEXPRG"].map({1.0: True, 2.0: False}).astype("boolean"),
        "sdmvpsu": d["SDMVPSU"].astype("Int64"), "sdmvstra": d["SDMVSTRA"].astype("Int64"),
        "wtint2yr": d["WTINT2YR"], "wtmec2yr": d["WTMEC2YR"], "wtmec4yr_pooled": d["WTMEC2YR"] / 2.0,
    })
    return out


def build_labs() -> pd.DataFrame:
    rows = []
    for role, vars_ in LAB_ROLE_VARS.items():
        for cyc, f in FILES[role].items():
            if not f:
                continue
            cb = codebook(f)
            lay = N.xport_layout(xpt(f))
            want = [HDL_VARS[cyc] if v == "HDL" else v for v in vars_]
            have = [v for v in want if v in lay["by_name"]]
            df = N.read_xport(xpt(f), ["SEQN"] + have)
            long = df.melt(id_vars="SEQN", var_name="lab_variable", value_name="value").dropna(subset=["value"])
            long["component_file"] = f
            long["component_role"] = role
            long["cycle_code"] = cyc
            long["lab_name"] = long["lab_variable"].map(lambda v: re.sub(r"\s*\([^()]*\)\s*$", "", cb.get(v, {}).get("label", v)))
            long["unit"] = long["lab_variable"].map(lambda v: cb.get(v, {}).get("unit", ""))
            long["subsample"] = LAB_SUBSAMPLE.get(role, "full MEC examined sample (weight WTMEC2YR)")
            long["harmonized_name"] = long["lab_variable"].map(lambda v: LAB_WIDE.get("HDL" if v in HDL_VARS.values() else v, ""))
            rows.append(long)
    labs = pd.concat(rows, ignore_index=True)
    labs["SEQN"] = labs["SEQN"].astype("int64")
    labs.insert(0, "participant_id", pid(labs["SEQN"]).values)
    labs["cycle"] = labs["cycle_code"].map(lambda c: CYCLES[c]["label"])
    labs["record_key"] = labs["component_file"] + ":SEQN=" + labs["SEQN"].astype(str) + ":" + labs["lab_variable"]
    return labs


def build_medications() -> pd.DataFrame:
    rx = load("RXQ_RX", ["RXDUSE", "RXDDRUG", "RXDCOUNT"])
    rx = rx.rename(columns={"RXDDRUG": "rx_generic_name", "RXDCOUNT": "rx_count_reported", "RXDUSE": "rx_use_code"})
    rx.insert(0, "participant_id", pid(rx["SEQN"]).values)
    rx["cycle"] = rx["cycle_code"].map(lambda c: CYCLES[c]["label"])
    rx["record_key"] = "RXQ_RX_" + rx["cycle_code"] + ":SEQN=" + rx["SEQN"].astype(str) + ":row=" + rx.index.astype(str)
    return rx


def build_mortality() -> pd.DataFrame:
    frames = []
    for cyc, fn in MORT_FILES.items():
        colspecs = [(0, 6), (14, 15), (15, 16), (16, 19), (19, 20), (20, 21), (42, 45), (45, 48)]
        names = ["SEQN", "mort_eligstat", "mort_status", "mort_ucod_leading_code", "mort_diabetes_mcod",
                 "mort_hyperten_mcod", "mort_permth_int", "mort_permth_exm"]
        df = pd.read_fwf(raw_dir(SOURCE_ID) / fn, colspecs=colspecs, names=names, na_values=["."], dtype=str)
        for c in names:
            df[c] = pd.to_numeric(df[c].str.strip(), errors="coerce")
        df["cycle_code"] = cyc
        frames.append(df)
    m = pd.concat(frames, ignore_index=True)
    m["SEQN"] = m["SEQN"].astype("int64")
    m["mort_ucod_leading"] = m["mort_ucod_leading_code"].map(N.MORT_UCOD)
    m["mort_eligstat"] = m["mort_eligstat"].map({1: "eligible", 2: "under age 18, not released", 3: "ineligible"})
    m["mort_followup_end"] = N.MORT_FOLLOWUP_END
    m.insert(0, "participant_id", pid(m["SEQN"]).values)
    m["cycle"] = m["cycle_code"].map(lambda c: CYCLES[c]["label"])
    m["record_key"] = m["cycle_code"].map(MORT_FILES) + ":SEQN=" + m["SEQN"].astype(str)
    return m


def build_clinical(parts: pd.DataFrame, labs: pd.DataFrame, meds: pd.DataFrame, mort: pd.DataFrame) -> pd.DataFrame:
    base = parts[["participant_id", "SEQN", "cycle", "cycle_code", "sex", "age_years", "race_ethnicity",
                  "education_adult20plus", "income_poverty_ratio", "mec_examined", "sdmvpsu", "sdmvstra",
                  "wtmec2yr"]].copy()
    seqn = base["SEQN"]
    wanted: dict[str, set] = {}
    for c, v, _ in YN_ITEMS:
        wanted.setdefault(c, set()).add(v)
    for c, v, *_ in RANGE_ITEMS:
        wanted.setdefault(c, set()).add(v)
    wanted.setdefault("BPX", set()).update({"BPXPLS", "BPXPULS"} | {f"BPXSY{i}" for i in range(1, 5)}
                                           | {f"BPXDI{i}" for i in range(1, 5)})
    wanted.setdefault("DIQ", set()).add("DIQ010")
    wanted.setdefault("MCQ", set()).update({"MCQ190", "MCQ191", "MCQ195"})
    wanted.setdefault("SMQ", set()).update({"SMQ020", "SMQ040"})
    wanted.setdefault("HUQ", set()).update({"HUQ050"})
    wanted.setdefault("PFQ", set()).update({"PFQ049", "PFQ051", "PFQ054", "PFQ057", "PFQ059"})
    wanted.setdefault("DPQ", set()).update(N.PHQ9_ITEMS)
    loaded = {}
    for comp, vs in wanted.items():
        if comp not in FILES:
            continue
        df = load(comp, sorted(vs))
        if comp == "RXQ_RX":
            df = df.groupby("SEQN", as_index=False).agg(RXDUSE=("RXDUSE", "first"))
        loaded[comp] = df.drop_duplicates("SEQN").set_index("SEQN")
    feats = {"SEQN": seqn}
    for comp, var, name in YN_ITEMS:
        if comp in loaded:
            feats[name] = N.yes_no(seqn.map(loaded[comp][var]))
    for comp, var, name, lo, hi in RANGE_ITEMS:
        if comp in loaded:
            feats[name] = N.in_range(seqn.map(loaded[comp][var]), lo, hi)
    bpx = loaded["BPX"]
    sy = bpx[[f"BPXSY{i}" for i in range(1, 5)]]
    di = bpx[[f"BPXDI{i}" for i in range(1, 5)]].where(lambda x: x > 0)
    feats["sbp_mean_mmhg"] = seqn.map(sy.mean(axis=1))
    feats["dbp_mean_mmhg"] = seqn.map(di.mean(axis=1))
    feats["exam_pulse_60s_bpm"] = seqn.map(bpx["BPXPLS"].where(bpx["BPXPLS"] > 0))
    pfq = loaded["PFQ"][["PFQ049", "PFQ051", "PFQ054", "PFQ057", "PFQ059"]].apply(N.yes_no)
    feats["any_functional_limitation_pfq"] = seqn.map(pfq.max(axis=1, skipna=True))
    diq = loaded["DIQ"]["DIQ010"]
    feats["dx_diabetes"] = seqn.map(diq.map({1.0: 1.0, 2.0: 0.0, 3.0: 0.0}))
    # arthritis type: 2003-2006 asks MCQ190 (1 rheumatoid, 2 osteoarthritis, 3 other) -- no psoriatic category
    mcq = loaded["MCQ"]
    at = mcq["MCQ190"].map({1.0: "rheumatoid arthritis", 2.0: "osteoarthritis", 3.0: "other"})
    if "MCQ191" in mcq:
        at = at.fillna(mcq["MCQ191"].map({1.0: "rheumatoid arthritis", 2.0: "osteoarthritis", 3.0: "psoriatic arthritis",
                                          4.0: "other"}))
    feats["arthritis_type"] = seqn.map(at)
    smq = loaded["SMQ"]
    status = np.select([smq["SMQ020"].eq(2.0), smq["SMQ020"].eq(1.0) & smq["SMQ040"].isin([1.0, 2.0]),
                        smq["SMQ020"].eq(1.0) & smq["SMQ040"].eq(3.0)], ["never", "current", "former"], default="")
    feats["smoking_status"] = seqn.map(pd.Series(status, index=smq.index).replace("", np.nan))
    if "DPQ" in loaded:
        phq = N._phq9(loaded["DPQ"].reset_index()).set_index("SEQN")
        for c in phq.columns:
            feats[c] = seqn.map(phq[c])
    mc = meds.dropna(subset=["rx_generic_name"]).groupby("SEQN").agg(rx_n_distinct_generic=("rx_generic_name", "nunique"))
    feats["rx_n_distinct_generic"] = seqn.map(mc["rx_n_distinct_generic"])
    feats["rx_n_distinct_generic"] = feats["rx_n_distinct_generic"].where(~feats["rx_any_past_month"].eq(0.0), 0)
    lw = labs.assign(key=labs["harmonized_name"]).query("key != ''")
    lw = lw.pivot_table(index="SEQN", columns="key", values="value", aggfunc="first")
    for name in LAB_WIDE.values():
        feats[f"lab_{name}"] = seqn.map(lw[name]) if name in lw.columns else pd.Series(np.nan, index=seqn.index)
    mm = mort.set_index("SEQN")
    for c in ["mort_eligstat", "mort_status", "mort_ucod_leading", "mort_permth_exm"]:
        feats[c] = seqn.map(mm[c])
    wide = base.merge(pd.DataFrame(feats), on="SEQN", how="left", validate="one_to_one")
    for c in ["dx_celiac_disease", "dx_copd", "told_sleep_disorder", "told_doctor_trouble_sleeping"]:
        if c not in wide.columns:
            wide[c] = np.nan
    wide["record_key"] = "SEQN=" + wide["SEQN"].astype(str)
    return wide


# --------------------------------------------------------------------------------------------- audit + registry


def _pct(n, d):
    return f"{n} ({100 * n / d:.1f}%)" if d else "0"


def write_audit(t: dict, wear: dict | None) -> Path:
    man = load_manifest(SOURCE_ID)["files"]
    P, C = t["participants"], t["clinical"]
    L = [f"# DATA AUDIT — NHANES 2003-2004 (C) and 2005-2006 (D): hip-accelerometry replication cohort", "",
         f"_Generated by `measure_it.ingestion.nhanes_2003_2006` (+ `measure_it.wearables.nhanes0306_features`) on "
         f"{utc_now_iso()}; every count below is computed from the retrieved files._", "",
         "| Field | Value |", "|---|---|",
         f"| source_id | {SOURCE_ID} |",
         "| Source | NHANES continuous survey cycles 2003-2004 (C) and 2005-2006 (D): DEMO, BPX, BMX, MCQ, BPQ, CDQ, DIQ, "
         "PFQ, HSQ, HUQ, SMQ, RXQ_RX, SLQ (D only), DPQ (D only), CBC/L25, BIOPRO/L40, GHB/L10, TCHOL+HDL/L13, VID, "
         "B12/L06NB, CRP/L11, ferritin (FERTIN_D / L06TFR_C), fasting GLU/L10AM and TRIGLY/L13AM, PAXRAW (hip "
         "accelerometer minute file); NCHS public-use Linked Mortality Files (2019) |",
         "| Publishing organization | CDC National Center for Health Statistics (NCHS) |",
         f"| Retrieval date (UTC) | {min(v['retrieved_at'] for v in man.values())} to {max(v['retrieved_at'] for v in man.values())} |",
         "| Source version / release | 2003-2004 and 2005-2006 public-use releases (per-file Last-Modified in MANIFEST.json; "
         "PAXRAW_C/D.zip last-modified 2014-12-31); LMF follow-up to 2019-12-31 |",
         f"| License / access conditions | {N.ACCESS_CONDITIONS} |",
         "| Unit of observation | participant (SEQN); accelerometry per participant-minute |",
         f"| Sample size (actual) | {len(P)} participants ({(P.cycle_code == 'C').sum()} C + {(P.cycle_code == 'D').sum()} D); "
         f"{int(P.mec_examined.sum())} MEC examined |",
         "| Geography | none released or inferred |", "| Person-level? | yes |", "| Geographic? | no |", "| Omics? | no |",
         "| Wearable? | yes — hip-worn ActiGraph AM-7164 uniaxial counts (+ steps in 2005-2006); NOT the 2011-2014 wrist MIMS |",
         "| True participant linkage across modalities? | yes, within NHANES on SEQN only; participant_id = nhanes0306:<SEQN> |",
         "", "## Files retrieved", "", "| file | bytes | sha256 (16) | Last-Modified | URL |", "|---|---:|---|---|---|"]
    for k in sorted(man):
        v = man[k]
        if k.startswith("docs/"):
            continue
        L.append(f"| {k} | {v['bytes']:,} | {v['sha256'][:16]} | {v.get('last_modified') or ''} | {v['url']} |")
    L += ["", f"Documentation pages saved under `docs/` ({sum(1 for k in man if k.startswith('docs/'))} files incl. "
          "listing pages).", "", "Requested but not published / not comparable:", ""]
    L += [f"* `{k}` — {v}" for k, v in NOT_PUBLISHED.items()]
    L += ["", "## Sample sizes and missingness (measured)", "", "| count | 2003-2004 (C) | 2005-2006 (D) | total |",
          "|---|---:|---:|---:|"]

    def row(label, s):
        c_, d_ = int(s[P.cycle_code == "C"].sum()), int(s[P.cycle_code == "D"].sum())
        L.append(f"| {label} | {c_} | {d_} | {c_ + d_} |")
    row("participants in DEMO", pd.Series(True, index=P.index))
    row("MEC examined", P["mec_examined"])
    row("aged >= 18", P["age_years"] >= 18)
    if wear:
        for k, v in wear.get("flow", {}).items():
            L.append(f"| {k} | {v.get('C', '')} | {v.get('D', '')} | {v.get('total', '')} |")
    ad = C[C["age_years"] >= 18]
    L += ["", "Missingness among adults (>= 18) for the analysis variables:", "",
          "| variable | n non-missing | missing (n, %) |", "|---|---:|---:|"]
    for c in [c for c in C.columns if c.startswith(("lab_", "dx_", "phq", "any_functional", "self_rated", "mort_",
                                                     "sbp", "dbp", "bmi", "waist", "exam_pulse", "rx_", "told_",
                                                     "smoking"))]:
        n = int(ad[c].notna().sum())
        L.append(f"| {c} | {n} | {_pct(len(ad) - n, len(ad))} |")
    if wear:
        L += ["", "## Hip accelerometry (PAXRAW) processing", ""] + wear.get("audit_lines", [])
    L += ["", "## Linkage strategy", "",
          "* Joins on SEQN only; participant_id = `nhanes0306:<SEQN>`. SEQN ranges: C 21005-31126, D 31127-41474 "
          "(measured in the tables); no SEQN is shared with 2011-2014 (62161-83731), and the datasets are kept "
          "separate anyway (different id namespace, different tables).",
          "* NOT joinable to any geographic, facility or other person-level dataset.", "",
          "## Limitations and caveats", "",
          "* Hip counts (ActiGraph AM-7164, proprietary counts/min, one axis) are a different device, placement and "
          "unit from the 2011-2014 wrist MIMS; no feature is treated as the same quantity. Any cross-cohort "
          "comparison is of conclusions, or of within-cohort standardised/percentile values labelled as such.",
          "* The device was worn during waking hours only (removed for sleep and water), so rest-activity metrics "
          "(IS/IV/M10/L5/RA) describe the waking-wear profile, not 24-h rhythms, and no sleep proxy is derived.",
          "* PHQ-9 (fatigue item, depression screen, ME/CFS-like proxy) exists only for 2005-2006; COPD and celiac "
          "questions do not exist in these cycles; arthritis type uses MCQ190 (no psoriatic category).",
          "* Age is top-coded at 85 in these cycles (80 in 2011-2014).",
          "* No heart rate/HRV stream; BPXPLS is an exam pulse.", "",
          "## Processed outputs", ""]
    from ..store import processed_path
    for name in ["participants__nhanes0306", "participant_clinical_features__nhanes0306", "participant_labs__nhanes0306",
                 "participant_medications__nhanes0306", "participant_mortality__nhanes0306",
                 "participant_wearable_features__nhanes0306", "participant_wearable_daily__nhanes0306"]:
        meta = processed_path(name).with_suffix(".meta.json")
        rows_ = json.loads(meta.read_text())["rows"] if meta.exists() else "not written"
        L.append(f"* `data/processed/{name}.parquet` — {rows_} rows")
    L += ["", "## Reproduce", "", "`uv run python -m measure_it.ingestion.nhanes_2003_2006` then "
          "`uv run python -m measure_it.wearables.nhanes0306_features`; tests `uv run pytest tests/test_nhanes_2003_2006.py`."]
    p = raw_dir(SOURCE_ID) / "DATA_AUDIT.md"
    p.write_text("\n".join(L) + "\n")
    return p


def write_registry(t: dict, wear: dict | None) -> Path:
    man = load_manifest(SOURCE_ID)["files"]
    P = t["participants"]
    entry = {
        "source_id": SOURCE_ID, "name": "NHANES 2003-2004 (C) and 2005-2006 (D) — hip accelerometry cohort",
        "publisher": "CDC National Center for Health Statistics (NCHS)", "landing_url": N.LANDING_URL,
        "access_urls": ["https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2003/DataFiles/",
                        "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2005/DataFiles/",
                        "https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality/"],
        "license": "U.S. government public-use data; NCHS Data User Agreement applies",
        "access_conditions": N.ACCESS_CONDITIONS, "retrieved_at": max(v["retrieved_at"] for v in man.values()),
        "source_version": "NHANES 2003-2004 and 2005-2006 public releases; PAXRAW_C/D.zip (2014-12-31); public-use LMF 2019",
        "update_date": N._latest_last_modified(man), "data_layer": "person",
        "unit_of_observation": "participant (SEQN); participant-minute for hip accelerometry",
        "sample_size": {"participants": int(len(P)), "per_cycle": {CYCLES[c]["label"]: int((P.cycle_code == c).sum())
                                                                   for c in CYCLES},
                        "mec_examined": int(P.mec_examined.sum()), **(wear or {}).get("registry_counts", {})},
        "geographic_resolution": "none", "person_level": True, "geographic": False, "omics": False, "wearable": True,
        "participant_linkage": "SEQN joins NHANES 2003-2006 components and the NCHS linked mortality file; "
                               "participant_id namespace nhanes0306:",
        "true_participant_linkage_across_modalities": True,
        "status": "ingested" if wear else "partial",
        "processed_outputs": ["participants__nhanes0306", "participant_clinical_features__nhanes0306",
                              "participant_labs__nhanes0306", "participant_medications__nhanes0306",
                              "participant_mortality__nhanes0306", "participant_wearable_features__nhanes0306",
                              "participant_wearable_daily__nhanes0306"],
        "audit": f"data/raw/{SOURCE_ID}/DATA_AUDIT.md",
        "ingestion_module": "measure_it.ingestion.nhanes_2003_2006 (+ measure_it.wearables.nhanes0306_features)",
        "limitations": ["Hip ActiGraph counts: different device/placement/unit from 2011-2014 wrist MIMS; not pooled",
                        "Waking-hours wear protocol: no sleep proxy; rest-activity metrics describe waking wear only",
                        "PHQ-9 (fatigue, depression, ME/CFS-like proxy) only 2005-2006",
                        "No heart rate/HRV stream", "No geography in public files"],
    }
    if not wear:
        entry["status_reason"] = "PAXRAW features not yet built"
    return write_registry_entry(entry)


def build_tables(log=print) -> dict:
    parts = build_participants()
    labs = build_labs()
    meds = build_medications()
    mort = build_mortality()
    clin = build_clinical(parts, labs, meds, mort)
    parts["record_key"] = "DEMO_" + parts["cycle_code"] + ":SEQN=" + parts["SEQN"].astype(str)
    write_person(parts, "participants__nhanes0306", evidence_type="person_self_report",
                 notes="Household-interview demographics + survey design variables. Age top-coded at 85.",
                 description="NHANES 2003-2006 participants (hip accelerometry replication cohort)")
    write_person(clin, "participant_clinical_features__nhanes0306", evidence_type="person_derived_feature",
                 notes="One row per DEMO participant assembled on SEQN; same variable names as the 2011-2014 table "
                       "where the item exists (absent items are NaN). BPXPLS is an exam pulse, not a wearable HR.",
                 description="NHANES 2003-2006 clinical features (wide)")
    write_person(labs, "participant_labs__nhanes0306", evidence_type="person_lab_measurement",
                 notes="Laboratory values as released (incl. CRP, ferritin subgroup, fasting subsample).",
                 description="NHANES 2003-2006 laboratory results (long)")
    write_person(meds.drop(columns=[]), "participant_medications__nhanes0306", evidence_type="person_self_report",
                 notes="Prescription medications, past 30 days (no reason-for-use codes in these cycles).",
                 description="NHANES 2003-2006 prescription medications (long)")
    write_person(mort, "participant_mortality__nhanes0306", evidence_type="person_linked_death_record",
                 notes="NCHS public-use Linked Mortality File (NDI linkage), follow-up through 2019-12-31.",
                 description="NHANES 2003-2006 linked mortality (public-use LMF 2019)")
    log(f"participants={len(parts)} clinical={len(clin)} labs={len(labs)} meds={len(meds)} mortality={len(mort)}")
    return {"participants": parts, "clinical": clin, "labs": labs, "medications": meds, "mortality": mort}


def run(log=print) -> dict:
    fetch_all(log=log)
    t = build_tables(log=log)
    write_audit(t, None)
    write_registry(t, None)
    return t


def _cli() -> None:
    argparse.ArgumentParser(description=__doc__.splitlines()[0]).parse_args()
    run()


if __name__ == "__main__":
    _cli()
