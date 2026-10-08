"""NHANES 2011-2012 (G) and 2013-2014 (H): person-level clinical + accelerometry ingestion.

Source files (verified 2026-09-23 against the NHANES data listing pages):
  * component XPT files:  https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/{2011|2013}/DataFiles/<FILE>.xpt
    with the documentation page next to each file (<FILE>.htm).
  * minute-level accelerometry (PAXMIN_G/H, 8.1 / 9.4 GB): the NHANES data listing links these to
    https://ftp.cdc.gov/pub/NHANES/LargeDataFiles/<FILE>.xpt; the DataFiles/ path returns HTTP 404.
  * public-use linked mortality files: https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality/

Joins are made ONLY on SEQN. participant_id = "nhanes:<SEQN>". Survey design variables
(SDMVPSU, SDMVSTRA, WTINT2YR, WTMEC2YR) are carried; pooling the two cycles requires WTMEC2YR/2.

Outputs (data/processed):
  participants__nhanes, participant_clinical_features__nhanes, participant_conditions__nhanes,
  participant_labs__nhanes, participant_medications__nhanes, participant_mortality__nhanes
  (wearable tables are built by measure_it.wearables.nhanes_features).

Run: uv run python -m measure_it.ingestion.nhanes [--skip-large]
"""
from __future__ import annotations

import argparse
import json
import os
import re
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
import requests

from ..config import RAW, UNKNOWN, raw_dir, utc_now_iso
from ..download import _locked, download_file, load_manifest, manifest_path, sha256_file
from ..http import USER_AGENT, HTMLInsteadOfData, ensure_online
from ..provenance import add_provenance
from ..registry import write_registry_entry
from ..store import write_table

SOURCE_ID = "nhanes_2011_2014"
SOURCE_NAME = "NHANES 2011-2012 (G) and 2013-2014 (H), CDC NCHS"
PRODUCER = "measure_it.ingestion.nhanes"
DATASET_PREFIX = "nhanes"

CYCLES = {
    "G": {"year": 2011, "label": "2011-2012", "sddsrvyr": 7},
    "H": {"year": 2013, "label": "2013-2014", "sddsrvyr": 8},
}
PUBLIC_BASE = "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/{year}/DataFiles/{file}"
LARGE_BASE = "https://ftp.cdc.gov/pub/NHANES/LargeDataFiles/{file}"
MORT_BASE = "https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality/{file}"
MORT_FILES = {"G": "NHANES_2011_2012_MORT_2019_PUBLIC.dat", "H": "NHANES_2013_2014_MORT_2019_PUBLIC.dat"}
MORT_READIN = "R_ReadInProgramAllSurveys.R"

# Components requested in SPEC §A.1. `cycles` lists the cycles in which the file exists on the
# NHANES data listing pages (checked 2026-09-23). Absent files are recorded, never substituted.
COMPONENTS: dict[str, dict] = {
    "DEMO": {"cycles": "GH", "kind": "demographics"},
    "BPX": {"cycles": "GH", "kind": "exam"},
    "BMX": {"cycles": "GH", "kind": "exam"},
    "MCQ": {"cycles": "GH", "kind": "questionnaire"},
    "BPQ": {"cycles": "GH", "kind": "questionnaire"},
    "CDQ": {"cycles": "GH", "kind": "questionnaire"},
    "DIQ": {"cycles": "GH", "kind": "questionnaire"},
    "SLQ": {"cycles": "GH", "kind": "questionnaire"},
    "PFQ": {"cycles": "GH", "kind": "questionnaire"},
    "DLQ": {"cycles": "H", "kind": "questionnaire"},
    "HSQ": {"cycles": "GH", "kind": "questionnaire"},
    "HUQ": {"cycles": "GH", "kind": "questionnaire"},
    "DPQ": {"cycles": "GH", "kind": "questionnaire"},
    "PAQ": {"cycles": "GH", "kind": "questionnaire"},
    "SMQ": {"cycles": "GH", "kind": "questionnaire"},
    "HIQ": {"cycles": "GH", "kind": "questionnaire"},
    "INQ": {"cycles": "GH", "kind": "questionnaire"},
    "RXQ_RX": {"cycles": "GH", "kind": "questionnaire"},
    "CBC": {"cycles": "GH", "kind": "lab"},
    "BIOPRO": {"cycles": "GH", "kind": "lab"},
    "GHB": {"cycles": "GH", "kind": "lab"},
    "GLU": {"cycles": "GH", "kind": "lab"},
    "TCHOL": {"cycles": "GH", "kind": "lab"},
    "HDL": {"cycles": "GH", "kind": "lab"},
    "TRIGLY": {"cycles": "GH", "kind": "lab"},
    "VID": {"cycles": "GH", "kind": "lab"},
    "VITB12": {"cycles": "GH", "kind": "lab"},
    "THYROD": {"cycles": "G", "kind": "lab"},
    "PAXHD": {"cycles": "GH", "kind": "accelerometry"},
    "PAXDAY": {"cycles": "GH", "kind": "accelerometry"},
    "PAXHR": {"cycles": "GH", "kind": "accelerometry"},
    "PAXMIN": {"cycles": "GH", "kind": "accelerometry_large"},
}
# Requested but not published for these cycles (probed 2026-09-23; see DATA_AUDIT.md).
NOT_PUBLISHED = {
    "FERTIN_G": "not on the 2011 laboratory data listing; DataFiles URL returns HTTP 404",
    "FERTIN_H": "not on the 2013 laboratory data listing; DataFiles URL returns HTTP 404",
    "HSCRP_G": "no CRP/hs-CRP file on the 2011 laboratory data listing; URL returns HTTP 404",
    "HSCRP_H": "no CRP/hs-CRP file on the 2013 laboratory data listing; URL returns HTTP 404",
    "CRP_G": "no CRP file on the 2011 laboratory data listing; URL returns HTTP 404",
    "CRP_H": "no CRP file on the 2013 laboratory data listing; URL returns HTTP 404",
    "THYROD_H": "thyroid profile not on the 2013 laboratory data listing; URL returns HTTP 404",
    "DLQ_G": "disability questionnaire (DLQ) starts in 2013-2014; URL returns HTTP 404",
    "PAXMIN_G (DataFiles path)": "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2011/DataFiles/PAXMIN_G.xpt returns HTTP 404; "
                                 "the data listing links https://ftp.cdc.gov/pub/NHANES/LargeDataFiles/PAXMIN_G.xpt",
    "PAXMIN_H (DataFiles path)": "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2013/DataFiles/PAXMIN_H.xpt returns HTTP 404; "
                                 "the data listing links https://ftp.cdc.gov/pub/NHANES/LargeDataFiles/PAXMIN_H.xpt",
}


# --------------------------------------------------------------------------------------------
# SAS XPORT (v5) reader: vectorised with numpy so multi-GB minute files stream in chunks without
# creating one Python object per character cell (pandas.read_sas does that for char columns).
# Validated against pandas.read_sas on several NHANES files in tests/test_nhanes.py.
# --------------------------------------------------------------------------------------------

def xport_layout(path: str | Path) -> dict:
    """Parse the header of a SAS XPORT v5 file: variables, row length, data offset, n rows."""
    path = Path(path)
    with open(path, "rb") as fh:
        head = fh.read(1 << 16)
    if not head.startswith(b"HEADER RECORD*******LIBRARY HEADER RECORD"):
        raise ValueError(f"{path.name}: not a SAS XPORT v5 file")
    mem = head.find(b"HEADER RECORD*******MEMBER  HEADER RECORD")
    namestr_len = int(head[mem + 74:mem + 78])
    nh = head.find(b"HEADER RECORD*******NAMESTR HEADER RECORD")
    nvars = int(head[nh + 54:nh + 58])
    start = nh + 80
    variables = []
    for i in range(nvars):
        rec = head[start + i * namestr_len:start + (i + 1) * namestr_len]
        ntype = int.from_bytes(rec[0:2], "big")
        length = int.from_bytes(rec[4:6], "big")
        name = rec[8:16].decode("ascii").strip()
        label = rec[16:56].decode("latin-1").strip()
        pos = int.from_bytes(rec[84:88], "big")
        variables.append({"name": name, "label": label, "numeric": ntype == 1, "length": length, "pos": pos})
    obs = head.find(b"HEADER RECORD*******OBS     HEADER RECORD", start + nvars * namestr_len - 80)
    data_start = obs + 80
    row_len = max(v["pos"] + v["length"] for v in variables)
    size = path.stat().st_size
    n_rows = (size - data_start) // row_len
    # trailing blank padding (to a multiple of 80 bytes) can look like one or more short final rows
    if n_rows and row_len < 80:
        with open(path, "rb") as fh:
            while n_rows:
                fh.seek(data_start + (n_rows - 1) * row_len)
                if fh.read(row_len).strip(b" ") != b"":
                    break
                n_rows -= 1
    return {"variables": variables, "row_len": row_len, "data_start": data_start, "n_rows": int(n_rows),
            "by_name": {v["name"].upper(): v for v in variables}}


def _ibm_to_float(raw: np.ndarray) -> np.ndarray:
    """IBM/360 hexadecimal floats (big-endian, 2..8 bytes, rows x bytes uint8) -> float64.

    SAS missing values ('.', '_', 'A'-'Z' in the first byte, zeros elsewhere) become NaN.
    """
    n, width = raw.shape
    buf = np.zeros((n, 8), dtype=np.uint8)
    buf[:, :width] = raw
    u = buf.view(">u8").ravel().astype(np.uint64)
    first = buf[:, 0]
    mant = u & np.uint64(0x00FFFFFFFFFFFFFF)
    expo = ((u >> np.uint64(56)) & np.uint64(0x7F)).astype(np.int64)
    val = np.ldexp(mant.astype(np.float64), (4 * (expo - 64) - 56).astype(np.int32))
    val = np.where((first & 0x80) != 0, -val, val)
    missing = (mant == 0) & (first != 0)
    val[missing] = np.nan
    return val


def iter_xport(path: str | Path, columns: list[str] | None = None, chunk_rows: int = 2_000_000,
               char_as: str = "str"):
    """Yield DataFrames of `chunk_rows` rows. Numeric -> float64; char -> str.

    char_as="bytes" keeps fixed-width bytes; char_as="uint8" returns 1-character columns as raw uint8
    byte codes (fast path for the minute files) and longer ones as bytes.
    """
    lay = xport_layout(path)
    cols = [c.upper() for c in (columns or [v["name"] for v in lay["variables"]])]
    missing = [c for c in cols if c not in lay["by_name"]]
    if missing:
        raise KeyError(f"{Path(path).name}: no variables {missing}")
    specs = [lay["by_name"][c] for c in cols]
    row_len, n_rows = lay["row_len"], lay["n_rows"]
    with open(path, "rb") as fh:
        fh.seek(lay["data_start"])
        done = 0
        while done < n_rows:
            k = min(chunk_rows, n_rows - done)
            block = np.frombuffer(fh.read(k * row_len), dtype=np.uint8)
            if block.size != k * row_len:
                raise IOError(f"{Path(path).name}: truncated at row {done}")
            rows = block.reshape(k, row_len)
            out = {}
            for v in specs:
                raw = rows[:, v["pos"]:v["pos"] + v["length"]]
                if v["numeric"]:
                    out[v["name"]] = _ibm_to_float(np.ascontiguousarray(raw))
                else:
                    if char_as == "uint8" and v["length"] == 1:
                        out[v["name"]] = raw[:, 0].copy()  # raw byte codes of 1-char columns
                        continue
                    s = np.ascontiguousarray(raw).view(f"S{v['length']}").ravel()
                    if char_as in ("bytes", "uint8"):
                        out[v["name"]] = s
                    else:
                        out[v["name"]] = np.char.strip(np.char.decode(s, "latin-1"))
            yield pd.DataFrame(out)
            done += k


def read_xport(path: str | Path, columns: list[str] | None = None) -> pd.DataFrame:
    """Read a whole (small) XPT file; empty character cells become NaN, numeric codes stay float."""
    parts = list(iter_xport(path, columns, chunk_rows=5_000_000))
    df = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=columns or [])
    for c in df.columns:
        if not pd.api.types.is_numeric_dtype(df[c]):
            df[c] = df[c].astype(object).where(df[c].astype(str) != "", np.nan)
    return df


def file_name(component: str, cycle: str) -> str:
    return f"{component}_{cycle}"


def xpt_url(component: str, cycle: str) -> str:
    fn = f"{file_name(component, cycle)}.xpt"
    if COMPONENTS[component]["kind"] == "accelerometry_large":
        return LARGE_BASE.format(file=fn)
    return PUBLIC_BASE.format(year=CYCLES[cycle]["year"], file=fn)


def doc_url(component: str, cycle: str) -> str:
    return PUBLIC_BASE.format(year=CYCLES[cycle]["year"], file=f"{file_name(component, cycle)}.htm")


def participant_id(seqn) -> pd.Series:
    s = pd.Series(seqn)
    return DATASET_PREFIX + ":" + s.astype("int64").astype(str)


# --------------------------------------------------------------------------------------------
# Downloads
# --------------------------------------------------------------------------------------------

def _record_manifest(rel: str, entry: dict) -> None:
    """Record a file downloaded outside download_file() in the shared MANIFEST.json (same lock)."""
    p = manifest_path(SOURCE_ID)
    with _locked(p):
        m = load_manifest(SOURCE_ID)
        m["files"][rel] = entry
        p.write_text(json.dumps(m, indent=2, sort_keys=True))


def download_large_ranged(url: str, filename: str, *, n_connections: int = 16,
                          chunk_bytes: int = 32 << 20, log=print) -> Path:
    """Segmented, resumable HTTP-range download for multi-GB files.

    ftp.cdc.gov serves ~100 kB/s per connection (measured 2026-09-23), so a single stream of
    the 8-9 GB PAXMIN files would take ~1 day. Chunks are written in place into <file>.part and
    finished chunk indices are kept in <file>.part.progress.json, so an interrupted run resumes.
    The finished file is recorded in MANIFEST.json with url, bytes, sha256 and retrieved_at.
    """
    dest = raw_dir(SOURCE_ID) / filename
    rel = filename
    manifest = load_manifest(SOURCE_ID)
    if dest.exists() and rel in manifest["files"]:
        return dest
    ensure_online(f"segmented download of {url}")
    hdrs = {"User-Agent": USER_AGENT}
    head = requests.head(url, headers=hdrs, timeout=60, allow_redirects=True)
    head.raise_for_status()
    if "text/html" in head.headers.get("Content-Type", ""):
        raise HTMLInsteadOfData(f"{url} is served as HTML")
    total = int(head.headers["Content-Length"])
    last_modified = head.headers.get("Last-Modified")
    part = dest.with_name(dest.name + ".part")
    prog_path = dest.with_name(dest.name + ".part.progress.json")
    done: set[int] = set()
    if part.exists() and prog_path.exists():
        prog = json.loads(prog_path.read_text())
        if prog.get("total") == total and prog.get("chunk_bytes") == chunk_bytes:
            done = set(prog["done"])
    if not part.exists() or not done:
        with open(part, "wb") as fh:
            fh.truncate(total)
    n_chunks = (total + chunk_bytes - 1) // chunk_bytes
    todo = [i for i in range(n_chunks) if i not in done]
    lock = threading.Lock()
    t0 = time.time()
    got = {"bytes": 0}

    def fetch(i: int) -> int:
        start = i * chunk_bytes
        end = min(total, start + chunk_bytes) - 1
        want = end - start + 1
        for attempt in range(8):
            try:
                buf = bytearray()
                with requests.get(url, headers={**hdrs, "Range": f"bytes={start}-{end}"},
                                  stream=True, timeout=120) as r:
                    if r.status_code != 206:
                        raise RuntimeError(f"range request answered HTTP {r.status_code}")
                    for block in r.iter_content(chunk_size=1 << 16):
                        buf.extend(block)
                if len(buf) != want:
                    raise RuntimeError(f"chunk {i}: got {len(buf)} of {want} bytes")
                fd = os.open(part, os.O_WRONLY)
                try:
                    os.pwrite(fd, bytes(buf), start)
                finally:
                    os.close(fd)
                with lock:
                    done.add(i)
                    got["bytes"] += want
                    prog_path.write_text(json.dumps({"url": url, "total": total, "chunk_bytes": chunk_bytes,
                                                     "done": sorted(done)}))
                return i
            except Exception as exc:  # network hiccup: back off and retry this chunk
                log(f"  chunk {i} attempt {attempt + 1} failed: {exc}")
                time.sleep(min(120, 5 * 2 ** attempt))
        raise RuntimeError(f"chunk {i} of {url} failed after retries")

    log(f"{filename}: {total} bytes, {n_chunks} chunks, {len(todo)} to fetch, {n_connections} connections")
    with ThreadPoolExecutor(max_workers=n_connections) as ex:
        futs = [ex.submit(fetch, i) for i in todo]
        for k, f in enumerate(as_completed(futs), 1):
            f.result()
            if k % 10 == 0 or k == len(todo):
                el = time.time() - t0
                log(f"  {filename}: {len(done)}/{n_chunks} chunks, {got['bytes'] / max(el, 1) / 1e6:.2f} MB/s")
    if len(done) != n_chunks or part.stat().st_size != total:
        raise RuntimeError(f"{filename}: incomplete download")
    os.replace(part, dest)
    prog_path.unlink(missing_ok=True)
    _record_manifest(rel, {
        "url": url, "final_url": url, "bytes": total, "sha256": sha256_file(dest),
        "retrieved_at": utc_now_iso(), "last_modified": last_modified, "etag": head.headers.get("ETag"),
        "content_type": head.headers.get("Content-Type"),
        "download_method": f"parallel HTTP range requests ({n_connections} connections, {chunk_bytes} B chunks)",
    })
    return dest


def fetch_component(component: str, cycle: str) -> Path:
    fn = f"{file_name(component, cycle)}.xpt"
    if COMPONENTS[component]["kind"] == "accelerometry_large":
        return download_large_ranged(xpt_url(component, cycle), fn)
    return download_file(xpt_url(component, cycle), SOURCE_ID, fn, min_bytes=1000)


def fetch_doc(component: str, cycle: str) -> Path:
    """Save the component's documentation page (HTML is expected here) under docs/."""
    (raw_dir(SOURCE_ID) / "docs").mkdir(exist_ok=True)
    return download_file(doc_url(component, cycle), SOURCE_ID, f"docs/{file_name(component, cycle)}.htm",
                         allow_html=True, min_bytes=5000)


def fetch_all(include_large: bool = True, log=print) -> dict[str, Path]:
    out: dict[str, Path] = {}
    for comp, meta in COMPONENTS.items():
        for cyc in meta["cycles"]:
            if meta["kind"] == "accelerometry_large" and not include_large:
                continue
            out[file_name(comp, cyc)] = fetch_component(comp, cyc)
            fetch_doc(comp, cyc)
    for cyc, fn in MORT_FILES.items():
        out[fn] = download_file(MORT_BASE.format(file=fn), SOURCE_ID, fn, min_bytes=1000)
    download_file(MORT_BASE.format(file=MORT_READIN), SOURCE_ID, MORT_READIN, min_bytes=100)
    return out




# --------------------------------------------------------------------------------------------
# Codebooks (documentation pages saved under docs/)
# --------------------------------------------------------------------------------------------

def _clean_html(s: str) -> str:
    import html as _html
    return re.sub(r"\s+", " ", _html.unescape(re.sub(r"<[^>]+>", "", s))).strip()


def parse_codebook(component: str, cycle: str, source_id: str = SOURCE_ID) -> dict[str, dict]:
    """Variable -> {label, unit, target, codes} parsed from the component's .htm documentation."""
    path = raw_dir(source_id) / "docs" / f"{file_name(component, cycle)}.htm"
    t = path.read_text(encoding="utf-8", errors="replace")
    out: dict[str, dict] = {}
    for blk in re.split(r'<div class="pagebreak">', t)[1:]:
        m = re.search(r'<h3 class="vartitle"[^>]*id="([^"]+)"[^>]*>(.*?)</h3>', blk, re.S)
        if not m:
            continue
        title = _clean_html(m.group(2))
        label = title.split(" - ", 1)[1] if " - " in title else title
        tgt = re.findall(r"Target: </dt>\s*<dd[^>]*>(.*?)</dd>", blk, re.S)
        rows = re.findall(r'<tr>\s*<td scope="row" class="values">(.*?)</td>\s*<td class="values">(.*?)</td>'
                          r'(?:\s*<td class="values" align="right">(.*?)</td>)?', blk, re.S)
        unit = re.search(r"\(([^()]*)\)\s*$", label)
        out[m.group(1).upper()] = {
            "label": label, "unit": unit.group(1).strip() if unit else "",
            "target": "; ".join(_clean_html(x) for x in tgt),
            "codes": {_clean_html(a): _clean_html(b) for a, b, _ in rows},
            "counts": {_clean_html(a): int(_clean_html(c)) for a, _, c in rows if _clean_html(c).isdigit()},
        }
    return out


# --------------------------------------------------------------------------------------------
# Loading and recoding helpers
# --------------------------------------------------------------------------------------------

def xpt_path(component: str, cycle: str) -> Path:
    return raw_dir(SOURCE_ID) / f"{file_name(component, cycle)}.xpt"


def load_component(component: str, columns: list[str] | None = None) -> pd.DataFrame:
    """Stack a component across the cycles in which it exists; adds cycle_code. Missing columns -> NaN."""
    frames = []
    for cyc in COMPONENTS[component]["cycles"]:
        path = xpt_path(component, cyc)
        lay = xport_layout(path)
        have = [c for c in (columns or [v["name"] for v in lay["variables"]]) if c.upper() in lay["by_name"]]
        df = read_xport(path, ["SEQN"] + [c for c in have if c.upper() != "SEQN"])
        for c in (columns or []):
            if c.upper() not in df.columns:
                df[c.upper()] = np.nan
        df["cycle_code"] = cyc
        frames.append(df)
    out = pd.concat(frames, ignore_index=True)
    out["SEQN"] = out["SEQN"].astype("int64")
    return out


def yes_no(s: pd.Series) -> pd.Series:
    """NHANES yes/no item: 1 -> 1.0, 2 -> 0.0, refused/don't know/missing -> NaN."""
    return s.map({1.0: 1.0, 2.0: 0.0}).astype("float64")


def in_range(s: pd.Series, lo: float, hi: float) -> pd.Series:
    """Keep values in [lo, hi]; NHANES refused/don't-know codes (7, 9, 77, 99, 777, ...) fall outside."""
    s = pd.to_numeric(s, errors="coerce")
    return s.where((s >= lo) & (s <= hi))


def cycle_meta(cycle: str) -> dict:
    m = load_manifest(SOURCE_ID)["files"].get(f"DEMO_{cycle}.xpt", {})
    return {
        "source_version": f"NHANES {CYCLES[cycle]['label']} public release (cycle {cycle}); "
                          f"DEMO_{cycle}.xpt last-modified {m.get('last_modified', UNKNOWN)}",
        "retrieved_at": m.get("retrieved_at", UNKNOWN),
    }


def _provenance_by_cycle(df: pd.DataFrame, *, evidence_type: str, record_id, notes: str,
                         evidence_level: str | pd.Series | None = None) -> pd.DataFrame:
    parts = []
    for cyc, sub in df.groupby("cycle_code", sort=True):
        meta = cycle_meta(cyc)
        lvl = evidence_level.loc[sub.index] if isinstance(evidence_level, pd.Series) else evidence_level
        parts.append(add_provenance(sub, data_layer="person", source_name=SOURCE_NAME,
                                    source_version=meta["source_version"], retrieved_at=meta["retrieved_at"],
                                    evidence_type=evidence_type, source_record_id=record_id,
                                    source_geographic_resolution="none", evidence_level=lvl,
                                    provenance_notes=notes))
    return pd.concat(parts, ignore_index=True)


# --------------------------------------------------------------------------------------------
# Participants (DEMO + PAXHD)
# --------------------------------------------------------------------------------------------

RACE3 = {1: "Mexican American", 2: "Other Hispanic", 3: "Non-Hispanic White", 4: "Non-Hispanic Black",
         6: "Non-Hispanic Asian", 7: "Other race including multi-racial"}
EDUC2 = {1: "Less than 9th grade", 2: "9-11th grade (no diploma)", 3: "High school graduate/GED",
         4: "Some college or AA degree", 5: "College graduate or above"}


def build_participants() -> pd.DataFrame:
    d = load_component("DEMO", ["SDDSRVYR", "RIDSTATR", "RIAGENDR", "RIDAGEYR", "RIDRETH3", "DMDEDUC2",
                                "INDFMPIR", "RIDEXMON", "RIDEXPRG", "SDMVPSU", "SDMVSTRA", "WTINT2YR", "WTMEC2YR"])
    pax = load_component("PAXHD", ["PAXSTS"])
    d = d.merge(pax[["SEQN", "PAXSTS"]], on="SEQN", how="left", validate="one_to_one")
    out = pd.DataFrame({
        "participant_id": participant_id(d["SEQN"]).values,
        "seqn": d["SEQN"].astype("int64"),
        "cycle": d["cycle_code"].map(lambda c: CYCLES[c]["label"]),
        "cycle_code": d["cycle_code"],
        "sddsrvyr": d["SDDSRVYR"].astype("Int64"),
        "interview_exam_status": d["RIDSTATR"].map({1.0: "interviewed only", 2.0: "interviewed and MEC examined"}),
        "mec_examined": d["RIDSTATR"].eq(2.0),
        "sex": d["RIAGENDR"].map({1.0: "male", 2.0: "female"}),
        "age_years": d["RIDAGEYR"],
        "age_topcoded_80": d["RIDAGEYR"].eq(80.0),
        "race_ethnicity": d["RIDRETH3"].map(RACE3),
        "education_adult20plus": d["DMDEDUC2"].map(EDUC2),
        "income_poverty_ratio": d["INDFMPIR"],
        "exam_period": d["RIDEXMON"].map({1.0: "Nov-Apr", 2.0: "May-Oct"}),
        "pregnant_at_exam": d["RIDEXPRG"].map({1.0: True, 2.0: False}),
        "sdmvpsu": d["SDMVPSU"].astype("Int64"),
        "sdmvstra": d["SDMVSTRA"].astype("Int64"),
        "wtint2yr": d["WTINT2YR"],
        "wtmec2yr": d["WTMEC2YR"],
        "wtmec4yr_pooled": d["WTMEC2YR"] / 2.0,
        "wtint4yr_pooled": d["WTINT2YR"] / 2.0,
        "pax_status": d["PAXSTS"].map({1.0: "PAM data available", 2.0: "no PAM data"}),
        "has_pam_data": d["PAXSTS"].eq(1.0),
    })
    out["pregnant_at_exam"] = out["pregnant_at_exam"].astype("boolean")
    return out


# --------------------------------------------------------------------------------------------
# Clinical features (one row per participant)
# --------------------------------------------------------------------------------------------

# (component, variable, output column, recode)
YN_ITEMS = [
    ("HUQ", "HUQ071", "hospitalized_overnight_past_year"), ("HUQ", "HUQ090", "mental_health_visit_past_year"),
    ("HIQ", "HIQ011", "health_insurance"),
    ("SLQ", "SLQ050", "told_doctor_trouble_sleeping"), ("SLQ", "SLQ060", "told_sleep_disorder"),
    ("PFQ", "PFQ049", "limitation_keeps_from_working"), ("PFQ", "PFQ051", "limited_amount_of_work"),
    ("PFQ", "PFQ054", "needs_special_equipment_to_walk"), ("PFQ", "PFQ057", "confusion_or_memory_problems"),
    ("PFQ", "PFQ059", "any_activity_limitation"), ("PFQ", "PFQ090", "needs_special_healthcare_equipment"),
    ("DLQ", "DLQ010", "serious_difficulty_hearing"), ("DLQ", "DLQ020", "serious_difficulty_seeing"),
    ("DLQ", "DLQ040", "serious_difficulty_concentrating"), ("DLQ", "DLQ050", "serious_difficulty_walking"),
    ("DLQ", "DLQ060", "difficulty_dressing_bathing"), ("DLQ", "DLQ080", "difficulty_errands_alone"),
    ("MCQ", "MCQ084", "confusion_memory_loss_worsening_past_year_60plus"),
    ("CDQ", "CDQ001", "chest_pain_ever_40plus"), ("CDQ", "CDQ010", "short_of_breath_stairs_inclines_40plus"),
    ("MCQ", "MCQ010", "dx_asthma"), ("MCQ", "MCQ160A", "dx_arthritis"), ("MCQ", "MCQ160N", "dx_gout"),
    ("MCQ", "MCQ160B", "dx_congestive_heart_failure"), ("MCQ", "MCQ160C", "dx_coronary_heart_disease"),
    ("MCQ", "MCQ160D", "dx_angina"), ("MCQ", "MCQ160E", "dx_heart_attack"), ("MCQ", "MCQ160F", "dx_stroke"),
    ("MCQ", "MCQ160G", "dx_emphysema"), ("MCQ", "MCQ160M", "dx_thyroid_problem"),
    ("MCQ", "MCQ160K", "dx_chronic_bronchitis"), ("MCQ", "MCQ160L", "dx_liver_condition"),
    ("MCQ", "MCQ160O", "dx_copd"), ("MCQ", "MCQ220", "dx_cancer"), ("MCQ", "MCQ082", "dx_celiac_disease"),
    ("MCQ", "MCQ070", "dx_psoriasis"), ("MCQ", "MCQ053", "anemia_treatment_past_3_months"),
    ("DIQ", "DIQ160", "dx_prediabetes"), ("DIQ", "DIQ050", "insulin_now"), ("DIQ", "DIQ070", "diabetes_pills_now"),
    ("BPQ", "BPQ020", "dx_hypertension"), ("BPQ", "BPQ080", "dx_high_cholesterol"),
    ("PAQ", "PAQ650", "vigorous_recreational_activity"), ("PAQ", "PAQ665", "moderate_recreational_activity"),
    ("RXQ_RX", "RXDUSE", "rx_any_past_month"),
]
# (component, variable, output column, lo, hi)
RANGE_ITEMS = [
    ("BMX", "BMXWT", "weight_kg", 0, 500), ("BMX", "BMXHT", "height_cm", 0, 300),
    ("BMX", "BMXBMI", "bmi", 0, 200), ("BMX", "BMXWAIST", "waist_cm", 0, 300),
    ("HSQ", "HSD010", "self_rated_health_mec_1excellent_5poor", 1, 5),
    ("HSQ", "HSQ470", "physically_unhealthy_days_30d", 0, 30),
    ("HSQ", "HSQ480", "mentally_unhealthy_days_30d", 0, 30),
    ("HSQ", "HSQ490", "inactive_days_due_to_health_30d", 0, 30),
    ("HSQ", "HSQ493", "pain_limited_activity_days_30d", 0, 30),
    ("HSQ", "HSQ496", "anxious_days_30d", 0, 30),
    ("HUQ", "HUQ010", "general_health_household_1excellent_5poor", 1, 5),
    ("HUQ", "HUD080", "n_overnight_hospital_stays_past_year", 1, 6),
    ("SLQ", "SLD010H", "sleep_hours_weekday_selfreport", 2, 12),
    ("MCQ", "MCQ380", "trouble_remembering_past_7d_60plus", 0, 4),
    ("PAQ", "PAD680", "sedentary_min_per_day_selfreport", 0, 1380),
    ("INQ", "INDFMMPI", "family_monthly_poverty_index", 0, 5),
    ("DPQ", "DPQ100", "phq_difficulty_caused_0to3", 0, 3),
]
PHQ9_ITEMS = [f"DPQ0{i}0" for i in range(1, 10)]
PFQ061_ITEMS = [f"PFQ061{c}" for c in "ABCDEFGHIJKLMNOPQRST"]

LAB_SPECS = {  # component -> NHANES variables kept in the long lab table (conventional units)
    "CBC": ["LBXWBCSI", "LBXLYPCT", "LBXMOPCT", "LBXNEPCT", "LBXEOPCT", "LBXBAPCT", "LBDLYMNO", "LBDMONO",
            "LBDNENO", "LBDEONO", "LBDBANO", "LBXRBCSI", "LBXHGB", "LBXHCT", "LBXMCVSI", "LBXMCHSI", "LBXMC",
            "LBXRDW", "LBXPLTSI", "LBXMPSI"],
    "BIOPRO": ["LBXSAL", "LBXSATSI", "LBXSASSI", "LBXSAPSI", "LBXSBU", "LBXSCA", "LBXSCK", "LBXSCH", "LBXSC3SI",
               "LBXSCR", "LBXSGTSI", "LBXSGL", "LBXSIR", "LBXSLDSI", "LBXSPH", "LBXSTB", "LBXSTP", "LBXSUA",
               "LBXSNASI", "LBXSKSI", "LBXSCLSI", "LBXSOSSI", "LBXSGB", "LBXSTR"],
    "GHB": ["LBXGH"],
    "GLU": ["LBXGLU", "LBXIN"],
    "TCHOL": ["LBXTC"],
    "HDL": ["LBDHDD"],
    "TRIGLY": ["LBXTR", "LBDLDL"],
    "VID": ["LBXVIDMS"],
    "VITB12": ["LBXB12", "LBDB12"],
    "THYROD": ["LBXTSH1", "LBXT4F", "LBXT3F", "LBXTT4", "LBXTT3", "LBXTPO", "LBXTGN", "LBXATG"],
}
LAB_SUBSAMPLE = {"GLU": "morning fasting subsample (weight WTSAF2YR)",
                 "TRIGLY": "morning fasting subsample (weight WTSAF2YR)",
                 "THYROD": "subsample A (weight WTSA2YR)"}
LAB_WIDE = {
    "LBXWBCSI": "wbc_1000_per_ul", "LBXLYPCT": "lymphocyte_pct", "LBXNEPCT": "neutrophil_pct",
    "LBXEOPCT": "eosinophil_pct", "LBXHGB": "hemoglobin_g_dl", "LBXHCT": "hematocrit_pct", "LBXMCVSI": "mcv_fl",
    "LBXRDW": "rdw_pct", "LBXPLTSI": "platelets_1000_per_ul", "LBXSAL": "albumin_g_dl", "LBXSATSI": "alt_u_l",
    "LBXSASSI": "ast_u_l", "LBXSAPSI": "alp_u_l", "LBXSBU": "bun_mg_dl", "LBXSCR": "creatinine_mg_dl",
    "LBXSCK": "ck_iu_l", "LBXSLDSI": "ldh_u_l", "LBXSUA": "uric_acid_mg_dl", "LBXSGL": "glucose_serum_mg_dl",
    "LBXSIR": "iron_ug_dl", "LBXSNASI": "sodium_mmol_l", "LBXSKSI": "potassium_mmol_l",
    "LBXSC3SI": "bicarbonate_mmol_l", "LBXSTP": "total_protein_g_dl", "LBXSGB": "globulin_g_dl",
    "LBXSCA": "calcium_mg_dl", "LBXSGTSI": "ggt_u_l", "LBXGH": "hba1c_pct", "LBXGLU": "glucose_fasting_mg_dl",
    "LBXTC": "total_cholesterol_mg_dl", "LBDHDD": "hdl_mg_dl", "LBXTR": "triglycerides_fasting_mg_dl",
    "LBDLDL": "ldl_fasting_mg_dl", "LBXVIDMS": "vitamin_d_25oh_nmol_l", "VITB12": "vitamin_b12_pg_ml",
    "LBXTSH1": "tsh_uiu_ml", "LBXT4F": "free_t4_ng_dl",
}


def build_labs() -> pd.DataFrame:
    rows = []
    for comp, vars_ in LAB_SPECS.items():
        for cyc in COMPONENTS[comp]["cycles"]:
            cb = parse_codebook(comp, cyc)
            lay = xport_layout(xpt_path(comp, cyc))
            have = [v for v in vars_ if v in lay["by_name"]]
            df = read_xport(xpt_path(comp, cyc), ["SEQN"] + have)
            long = df.melt(id_vars="SEQN", var_name="lab_variable", value_name="value").dropna(subset=["value"])
            long["component_file"] = file_name(comp, cyc)
            long["cycle_code"] = cyc
            long["lab_name"] = long["lab_variable"].map(lambda v: re.sub(r"\s*\([^()]*\)\s*$", "", cb.get(v, {}).get("label", v)))
            long["unit"] = long["lab_variable"].map(lambda v: cb.get(v, {}).get("unit", ""))
            long["subsample"] = LAB_SUBSAMPLE.get(comp, "full MEC examined sample (weight WTMEC2YR)")
            note = ""
            if comp == "VITB12" and cyc == "H":
                note = "LBDB12: 2013-2014 values adjusted by NCHS with Deming regression for reagent-lot QC shifts"
            long["lab_note"] = note
            rows.append(long)
    labs = pd.concat(rows, ignore_index=True)
    labs["SEQN"] = labs["SEQN"].astype("int64")
    labs.insert(0, "participant_id", participant_id(labs["SEQN"]).values)
    labs["cycle"] = labs["cycle_code"].map(lambda c: CYCLES[c]["label"])
    labs["harmonized_name"] = labs["lab_variable"].map(lambda v: LAB_WIDE.get("VITB12" if v in ("LBXB12", "LBDB12") else v, ""))
    labs["record_key"] = labs["component_file"] + ":SEQN=" + labs["SEQN"].astype(str) + ":" + labs["lab_variable"]
    return labs


def _phq9(dpq: pd.DataFrame) -> pd.DataFrame:
    items = pd.DataFrame({c: in_range(dpq[c], 0, 3) for c in PHQ9_ITEMS})
    n = items.notna().sum(axis=1)
    total = items.sum(axis=1).where(n == 9)
    return pd.DataFrame({"SEQN": dpq["SEQN"], "phq9_total": total, "phq9_n_items": n,
                         "phq9_ge10": (total >= 10).where(total.notna()).astype("float64"),
                         "phq9_item3_sleep_0to3": items["DPQ030"],
                         "phq9_item4_tired_little_energy_0to3": items["DPQ040"]})


def build_clinical(participants: pd.DataFrame, labs: pd.DataFrame, meds: pd.DataFrame,
                   mortality: pd.DataFrame | None) -> pd.DataFrame:
    base = participants[["participant_id", "seqn", "cycle", "cycle_code", "sex", "age_years", "age_topcoded_80",
                         "race_ethnicity", "education_adult20plus", "income_poverty_ratio", "mec_examined",
                         "sdmvpsu", "sdmvstra", "wtint2yr", "wtmec2yr", "has_pam_data"]].copy()
    base = base.rename(columns={"seqn": "SEQN"})
    comps = sorted({c for c, *_ in YN_ITEMS} | {c for c, *_ in RANGE_ITEMS} | {"BPX", "DIQ", "MCQ", "SMQ", "PFQ", "DPQ"})
    wanted: dict[str, set] = {c: set() for c in comps}
    for c, v, _ in YN_ITEMS:
        wanted[c].add(v)
    for c, v, *_ in RANGE_ITEMS:
        wanted[c].add(v)
    wanted["BPX"] |= {"BPXPLS", "BPXPULS"} | {f"BPXSY{i}" for i in range(1, 5)} | {f"BPXDI{i}" for i in range(1, 5)}
    wanted["DIQ"] |= {"DIQ010"}
    wanted["MCQ"] |= {"MCQ195"}
    wanted["SMQ"] |= {"SMQ020", "SMQ040"}
    wanted["HUQ"] |= {"HUQ050", "HUQ051"}
    wanted["PFQ"] |= set(PFQ061_ITEMS) | {"PFQ049", "PFQ051", "PFQ054", "PFQ057", "PFQ059"}
    wanted["DPQ"] |= set(PHQ9_ITEMS)
    wide = base
    loaded: dict[str, pd.DataFrame] = {}
    for comp in comps:
        df = load_component(comp, sorted(wanted[comp]))
        if comp == "RXQ_RX":  # one row per medication: participant-level items only
            df = df.groupby("SEQN", as_index=False).agg(RXDUSE=("RXDUSE", "first"))
        loaded[comp] = df
    seqn = base["SEQN"]
    feats: dict[str, pd.Series] = {"SEQN": seqn}
    for comp, var, name in YN_ITEMS:
        src = loaded[comp].set_index("SEQN")[var]
        feats[name] = yes_no(seqn.map(src))
    for comp, var, name, lo, hi in RANGE_ITEMS:
        src = loaded[comp].set_index("SEQN")[var]
        feats[name] = in_range(seqn.map(src), lo, hi)
    # blood pressure / exam pulse (BPXPLS is a 60-s radial pulse counted at the MEC exam; NOT a wearable HR)
    bpx = loaded["BPX"].set_index("SEQN")
    sy = bpx[[f"BPXSY{i}" for i in range(1, 5)]]
    di = bpx[[f"BPXDI{i}" for i in range(1, 5)]].where(lambda x: x > 0)  # diastolic 0 recorded -> missing
    feats["sbp_mean_mmhg"] = seqn.map(sy.mean(axis=1))
    feats["dbp_mean_mmhg"] = seqn.map(di.mean(axis=1))
    feats["n_sbp_readings"] = seqn.map(sy.notna().sum(axis=1)).fillna(0).astype(int)
    feats["exam_pulse_60s_bpm"] = seqn.map(bpx["BPXPLS"].where(bpx["BPXPLS"] > 0))
    feats["exam_pulse_irregular"] = seqn.map(bpx["BPXPULS"].map({1.0: 0.0, 2.0: 1.0}))
    # healthcare visits past year: HUQ050 (2011-12: 0 none,1,2-3,4-9,10-12,13+) vs HUQ051 (2013-14: finer bins)
    huq = loaded["HUQ"].set_index("SEQN")
    g_cat = in_range(huq["HUQ050"], 0, 5)
    h_cat = in_range(huq["HUQ051"], 0, 8).map({0: 0, 1: 1, 2: 2, 3: 3, 4: 3, 5: 3, 6: 4, 7: 5, 8: 5})
    feats["healthcare_visits_past_year_cat_0none_1one_2two3_3four9_4ten12_5thirteenplus"] = seqn.map(g_cat.fillna(h_cat))
    # any functional limitation: PFQ059 is only asked when PFQ049/051/054/057 are all 'no' (check item PFQ.058)
    pfq_any = loaded["PFQ"].set_index("SEQN")[["PFQ049", "PFQ051", "PFQ054", "PFQ057", "PFQ059"]].apply(yes_no)
    feats["any_functional_limitation_pfq"] = seqn.map(pfq_any.max(axis=1, skipna=True))
    # diabetes (DIQ010: 1 yes, 2 no, 3 borderline)
    diq = loaded["DIQ"].set_index("SEQN")["DIQ010"]
    feats["dx_diabetes"] = seqn.map(diq.map({1.0: 1.0, 2.0: 0.0, 3.0: 0.0}))
    feats["dx_borderline_diabetes"] = seqn.map(diq.map({1.0: 0.0, 2.0: 0.0, 3.0: 1.0}))
    mcq195 = loaded["MCQ"].set_index("SEQN")["MCQ195"]
    feats["arthritis_type"] = seqn.map(mcq195.map({1.0: "osteoarthritis", 2.0: "rheumatoid arthritis",
                                                              3.0: "psoriatic arthritis", 4.0: "other"}))
    smq = loaded["SMQ"].set_index("SEQN")
    status = np.select([smq["SMQ020"].eq(2.0), smq["SMQ020"].eq(1.0) & smq["SMQ040"].isin([1.0, 2.0]),
                        smq["SMQ020"].eq(1.0) & smq["SMQ040"].eq(3.0)], ["never", "current", "former"], default="")
    feats["smoking_status"] = seqn.map(pd.Series(status, index=smq.index).replace("", np.nan))
    # PFQ061 difficulty items: 1 none .. 4 unable; 5 'do not do this activity' and 7/9 -> NaN
    pfq = loaded["PFQ"].set_index("SEQN")
    pf = pd.DataFrame({c: in_range(pfq[c], 1, 4) for c in PFQ061_ITEMS})
    for c in PFQ061_ITEMS:
        feats[f"{c.lower()}_difficulty_1none_4unable"] = seqn.map(pf[c])
    feats["pfq061_n_items_answered"] = seqn.map(pf.notna().sum(axis=1))
    feats["pfq061_n_items_some_difficulty_or_worse"] = seqn.map((pf >= 2).sum(axis=1)).where(
        feats["pfq061_n_items_answered"] > 0)
    phq = _phq9(loaded["DPQ"]).set_index("SEQN")
    for c in phq.columns:
        feats[c] = seqn.map(phq[c])
    # medication counts
    mc = meds.groupby("SEQN").agg(rx_n_records=("rx_generic_name", "size"),
                                   rx_n_distinct_generic=("rx_generic_name", "nunique"),
                                   rx_count_reported=("rx_count_reported", "max"))
    for c in mc.columns:
        feats[c] = seqn.map(mc[c])
    no_rx = feats["rx_any_past_month"].eq(0.0)
    for c in ["rx_n_records", "rx_n_distinct_generic", "rx_count_reported"]:
        feats[c] = feats[c].where(~no_rx, 0)
    # wide labs
    lw = labs.assign(key=np.where(labs["lab_variable"].isin(["LBXB12", "LBDB12"]), "VITB12", labs["lab_variable"]))
    lw = lw[lw["key"].isin(LAB_WIDE)].pivot_table(index="SEQN", columns="key", values="value", aggfunc="first")
    for k, name in LAB_WIDE.items():
        feats[f"lab_{name}"] = seqn.map(lw[k]) if k in lw.columns else pd.Series(np.nan, index=seqn.index)
    if mortality is not None:
        mm = mortality.set_index("SEQN")
        for c in ["mort_eligstat", "mort_status", "mort_ucod_leading", "mort_permth_exm", "mort_permth_int"]:
            feats[c] = seqn.map(mm[c])
    feats = pd.DataFrame(feats)
    wide = base.merge(feats, on="SEQN", how="left", validate="one_to_one")
    wide["record_key"] = "SEQN=" + wide["SEQN"].astype(str)
    return wide.rename(columns={"SEQN": "seqn"})


# --------------------------------------------------------------------------------------------
# Conditions (long): self-reported diagnoses, questionnaire screens, rx reason-for-use codes
# --------------------------------------------------------------------------------------------

SELF_REPORT_CONDITIONS = [  # (component, variable, positive code, condition label)
    ("MCQ", "MCQ010", 1.0, "asthma"), ("MCQ", "MCQ160A", 1.0, "arthritis"),
    ("MCQ", "MCQ160N", 1.0, "gout"), ("MCQ", "MCQ160B", 1.0, "congestive heart failure"),
    ("MCQ", "MCQ160C", 1.0, "coronary heart disease"), ("MCQ", "MCQ160D", 1.0, "angina"),
    ("MCQ", "MCQ160E", 1.0, "heart attack"), ("MCQ", "MCQ160F", 1.0, "stroke"),
    ("MCQ", "MCQ160G", 1.0, "emphysema"), ("MCQ", "MCQ160M", 1.0, "thyroid problem"),
    ("MCQ", "MCQ160K", 1.0, "chronic bronchitis"), ("MCQ", "MCQ160L", 1.0, "liver condition"),
    ("MCQ", "MCQ160O", 1.0, "COPD"), ("MCQ", "MCQ220", 1.0, "cancer or malignancy"),
    ("MCQ", "MCQ082", 1.0, "celiac disease"), ("MCQ", "MCQ070", 1.0, "psoriasis"),
    ("MCQ", "MCQ195", 1.0, "osteoarthritis"), ("MCQ", "MCQ195", 2.0, "rheumatoid arthritis"),
    ("MCQ", "MCQ195", 3.0, "psoriatic arthritis"), ("MCQ", "MCQ195", 4.0, "other arthritis"),
    ("MCQ", "MCQ053", 1.0, "anemia treatment in past 3 months"),
    ("DIQ", "DIQ010", 1.0, "diabetes"), ("DIQ", "DIQ010", 3.0, "borderline diabetes"),
    ("DIQ", "DIQ160", 1.0, "prediabetes"), ("BPQ", "BPQ020", 1.0, "hypertension"),
    ("BPQ", "BPQ080", 1.0, "high cholesterol"),
    ("SLQ", "SLQ060", 1.0, "sleep disorder"),
    ("SLQ", "SLQ050", 1.0, "trouble sleeping (told a doctor)"),
]
RX_BASIS = "condition inferred from prescription reason-for-use code"
# ICD-10-CM prefixes for the invisible-illness conditions in configs/SPEC; NHANES appends 'P' for
# prevention use (e.g. G43.P 'Prevent migraine'), which the prefix match keeps. Every prefix was checked
# against the CMS ICD-10-CM table (ontology_icd10cm_codes) and the RXQ_RX_H codebook appendix. "chronic pain"
# is G89.2x (Chronic pain, NEC) and G89.4 (Chronic pain syndrome) only: the G89 category is "Pain, not
# elsewhere classified" and also holds acute pain (G89.1x, e.g. G89.18 'Other acute postprocedural pain',
# the most frequent G89 reason code in RXQ_RX_H) and G89.3 'Neoplasm related pain (acute) (chronic)'.
RX_TARGET_GROUPS = {
    "fibromyalgia": ["M79.7"], "migraine": ["G43"], "irritable bowel syndrome": ["K58"],
    "fatigue / malaise": ["R53"], "postviral fatigue syndrome": ["G93.3"],
    "orthostatic hypotension": ["I95.1"], "syncope": ["R55"], "autonomic nervous system disorder": ["G90"],
    "chronic pain": ["G89.2", "G89.4"], "myalgia": ["M79.1"], "endometriosis": ["N80"], "gastroparesis": ["K31.84"],
    "Ehlers-Danlos syndrome": ["Q79.6"], "Lyme disease": ["A69.2"], "tachycardia": ["R00.0"],
    "headache": ["R51"], "insomnia": ["G47.0", "F51.0"],
}


def rx_target_group(code: str) -> str:
    if not isinstance(code, str):
        return ""
    for grp, prefixes in RX_TARGET_GROUPS.items():
        if any(code.startswith(p) for p in prefixes):
            return grp
    return ""


def build_conditions(clinical: pd.DataFrame, meds: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for comp, var, code, label in SELF_REPORT_CONDITIONS:
        df = load_component(comp, [var])
        hit = df[df[var] == code]
        cb = parse_codebook(comp, hit["cycle_code"].iloc[0]) if len(hit) else {}
        rows.append(pd.DataFrame({
            "SEQN": hit["SEQN"], "cycle_code": hit["cycle_code"], "condition_label": label,
            "source_variable": f"{comp}.{var}", "value": f"{int(code)} = " + cb.get(var, {}).get("codes", {}).get(str(int(code)), ""),
            "icd10cm_code": "", "evidence_basis": "self-reported health-professional diagnosis"
            if var not in ("SLQ050", "MCQ053") else "self-reported item",
            "question": cb.get(var, {}).get("label", ""),
        }))
    # questionnaire-derived symptom screens (proxies; not diagnoses)
    c = clinical
    phq = c[c["phq9_ge10"] == 1.0]
    rows.append(pd.DataFrame({"SEQN": phq["seqn"], "cycle_code": phq["cycle_code"],
                              "condition_label": "depressive symptoms (PHQ-9 total >= 10)",
                              "source_variable": "DPQ.DPQ010-DPQ090", "value": phq["phq9_total"].map(lambda v: f"PHQ-9 total = {v:.0f}"),
                              "icd10cm_code": "", "evidence_basis": "questionnaire screen (symptom proxy, not a diagnosis)",
                              "question": "PHQ-9 depression screener"}))
    fat = c[c["phq9_item4_tired_little_energy_0to3"] >= 2]
    rows.append(pd.DataFrame({"SEQN": fat["seqn"], "cycle_code": fat["cycle_code"],
                              "condition_label": "fatigue symptom (tired / little energy more than half the days)",
                              "source_variable": "DPQ.DPQ040",
                              "value": fat["phq9_item4_tired_little_energy_0to3"].map({2.0: "2 = More than half the days", 3.0: "3 = Nearly every day"}),
                              "icd10cm_code": "", "evidence_basis": "self-reported symptom item (fatigue proxy, not a diagnosis)",
                              "question": "Feeling tired or having little energy (past 2 weeks)"}))
    sr = pd.concat(rows, ignore_index=True)
    sr["rx_target_group"] = ""
    # prescription reason-for-use ICD-10-CM codes (2013-2014 only)
    rx = []
    for i in (1, 2, 3):
        sub = meds[["SEQN", "cycle_code", f"reason_icd10cm_{i}", f"reason_description_{i}", "rx_generic_name"]].copy()
        sub.columns = ["SEQN", "cycle_code", "icd10cm_code", "desc", "drug"]
        sub["source_variable"] = f"RXQ_RX.RXDRSC{i}"
        rx.append(sub)
    rx = pd.concat(rx, ignore_index=True).dropna(subset=["icd10cm_code"])
    rx = rx[~rx["icd10cm_code"].isin(["55555", "77777", "99999"])]
    rx = (rx.groupby(["SEQN", "cycle_code", "icd10cm_code"], as_index=False)
            .agg(desc=("desc", "first"), drugs=("drug", lambda s: "; ".join(sorted(set(map(str, s))))),
                 source_variable=("source_variable", lambda s: ",".join(sorted(set(s))))))
    rx_rows = pd.DataFrame({
        "SEQN": rx["SEQN"], "cycle_code": rx["cycle_code"], "condition_label": rx["desc"],
        "source_variable": rx["source_variable"], "value": "prescribed: " + rx["drugs"],
        "icd10cm_code": rx["icd10cm_code"], "evidence_basis": RX_BASIS,
        "question": "Main reason for use of a prescription medication taken in the past 30 days (coded by NCHS)",
        "rx_target_group": rx["icd10cm_code"].map(rx_target_group),
    })
    out = pd.concat([sr, rx_rows], ignore_index=True)
    out["SEQN"] = out["SEQN"].astype("int64")
    out.insert(0, "participant_id", participant_id(out["SEQN"]).values)
    out["cycle"] = out["cycle_code"].map(lambda c: CYCLES[c]["label"])
    out["record_key"] = (out["source_variable"] + ":SEQN=" + out["SEQN"].astype(str) + ":"
                         + out["condition_label"].astype(str).str.slice(0, 60))
    return out


# --------------------------------------------------------------------------------------------
# Medications (long) and mortality
# --------------------------------------------------------------------------------------------

def build_medications() -> pd.DataFrame:
    cols = ["RXDUSE", "RXDDRUG", "RXDDRGID", "RXQSEEN", "RXDDAYS", "RXDCOUNT"] + \
           [f"RXDRSC{i}" for i in (1, 2, 3)] + [f"RXDRSD{i}" for i in (1, 2, 3)]
    rx = load_component("RXQ_RX", cols)
    rx = rx[rx["RXDUSE"].eq(1.0) & rx["RXDDRUG"].notna()].copy()
    special = {"55555": "unknown", "77777": "refused", "99999": "don't know"}
    out = pd.DataFrame({
        "SEQN": rx["SEQN"].astype("int64"), "cycle_code": rx["cycle_code"],
        "rx_generic_name": rx["RXDDRUG"].astype(str),
        "rx_name_status": rx["RXDDRUG"].astype(str).map(special).fillna("recorded"),
        # RXQSEEN (codebook): 1 Yes, 2 No, 3 'Only pharmacy print out seen' (2013-2014 only; container not seen)
        "rx_drug_code": rx["RXDDRGID"],
        "container_seen": rx["RXQSEEN"].map({1.0: True, 2.0: False, 3.0: False}).astype("boolean"),
        "container_seen_detail": rx["RXQSEEN"].map({1.0: "container seen", 2.0: "not seen",
                                                    3.0: "only pharmacy print out seen"}),
        "days_taken": in_range(rx["RXDDAYS"], 1, 50000), "rx_count_reported": rx["RXDCOUNT"],
    })
    for i in (1, 2, 3):
        out[f"reason_icd10cm_{i}"] = rx[f"RXDRSC{i}"].values
        out[f"reason_description_{i}"] = rx[f"RXDRSD{i}"].values
    out["reason_target_groups"] = [
        "; ".join(sorted({g for g in (rx_target_group(c) for c in codes) if g}))
        for codes in zip(out["reason_icd10cm_1"], out["reason_icd10cm_2"], out["reason_icd10cm_3"])]
    out = out.reset_index(drop=True)
    out.insert(0, "participant_id", participant_id(out["SEQN"]).values)
    out["cycle"] = out["cycle_code"].map(lambda c: CYCLES[c]["label"])
    out["rx_record_index"] = out.groupby("SEQN").cumcount() + 1
    out["record_key"] = "RXQ_RX_" + out["cycle_code"] + ":SEQN=" + out["SEQN"].astype(str) + ":rx" + out["rx_record_index"].astype(str)
    return out


MORT_UCOD = {1: "Diseases of heart", 2: "Malignant neoplasms", 3: "Chronic lower respiratory diseases",
             4: "Accidents (unintentional injuries)", 5: "Cerebrovascular diseases", 6: "Alzheimer's disease",
             7: "Diabetes mellitus", 8: "Influenza and pneumonia",
             9: "Nephritis, nephrotic syndrome and nephrosis", 10: "All other causes (residual)"}
MORT_FOLLOWUP_END = "2019-12-31"


def build_mortality() -> pd.DataFrame:
    """Public-use Linked Mortality File (NDI linkage through 2019-12-31), layout per NCHS R read-in program."""
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
    m["mort_ucod_leading"] = m["mort_ucod_leading_code"].map(MORT_UCOD)
    m["mort_eligstat"] = m["mort_eligstat"].map({1: "eligible", 2: "under age 18, not released", 3: "ineligible"})
    m["mort_followup_end"] = MORT_FOLLOWUP_END
    m.insert(0, "participant_id", participant_id(m["SEQN"]).values)
    m["cycle"] = m["cycle_code"].map(lambda c: CYCLES[c]["label"])
    m["record_key"] = m["cycle_code"].map(MORT_FILES) + ":SEQN=" + m["SEQN"].astype(str)
    return m


# --------------------------------------------------------------------------------------------
# Orchestration
# --------------------------------------------------------------------------------------------

def write_person_table(df: pd.DataFrame, name: str, *, evidence_type: str, notes: str,
                       description: str, evidence_level=None) -> Path:
    if "cycle_code" not in df.columns:
        raise ValueError("person tables need cycle_code for provenance")
    out = _provenance_by_cycle(df, evidence_type=evidence_type, record_id="record_key", notes=notes,
                               evidence_level=evidence_level)
    out = out.drop(columns=["record_key"]).rename(columns={"SEQN": "seqn"})
    return write_table(out, name, producer=PRODUCER, description=description)


def build_clinical_tables(log=print) -> dict[str, pd.DataFrame]:
    participants = build_participants()
    labs = build_labs()
    meds = build_medications()
    mort = build_mortality()
    clinical = build_clinical(participants, labs, meds, mort)
    conditions = build_conditions(clinical, meds)
    participants["record_key"] = "DEMO_" + participants["cycle_code"] + ":SEQN=" + participants["seqn"].astype(str)
    write_person_table(participants, "participants__nhanes", evidence_type="person_self_report",
                       notes="Household-interview demographics plus NHANES survey design variables "
                             "(SDMVPSU, SDMVSTRA, WTINT2YR, WTMEC2YR). Pool two cycles with WTMEC2YR/2. "
                             "No geography is released or inferred.",
                       description="NHANES 2011-2014 participants: demographics, survey design, cycle, PAM status")
    write_person_table(clinical, "participant_clinical_features__nhanes", evidence_type="person_derived_feature",
                       notes="One row per DEMO participant assembled on SEQN from DEMO, BPX, BMX, MCQ, BPQ, CDQ, DIQ, "
                             "SLQ, PFQ, DLQ (2013-14 only), HSQ, HUQ, DPQ, PAQ, SMQ, HIQ, INQ, RXQ_RX, labs and the "
                             "public-use linked mortality file. exam_pulse_60s_bpm is a 60-s radial pulse at the MEC "
                             "exam, not a wearable heart rate. Refused/don't know -> NaN.",
                       description="NHANES 2011-2014 clinical features (wide, one row per participant)")
    write_person_table(conditions, "participant_conditions__nhanes", evidence_type="person_self_report",
                       evidence_level=conditions["evidence_basis"],
                       notes="Affirmative answers only. evidence_basis distinguishes self-reported diagnoses, "
                             f"questionnaire symptom screens (proxies) and '{RX_BASIS}' (RXDRSC1-3, 2013-2014 only).",
                       description="NHANES 2011-2014 conditions (long)")
    write_person_table(labs, "participant_labs__nhanes", evidence_type="person_lab_measurement",
                       notes="Laboratory values as released (below-LOD values carry NCHS fill values LLOD/sqrt(2)); "
                             "fasting and subsample A analytes need their subsample weights.",
                       description="NHANES 2011-2014 laboratory results (long)")
    write_person_table(meds, "participant_medications__nhanes", evidence_type="person_self_report",
                       notes="Prescription medications taken in the past 30 days (household interview; container seen "
                             "flag). Reason-for-use ICD-10-CM codes exist only in 2013-2014 (RXQ_RX_H).",
                       description="NHANES 2011-2014 prescription medications (long)")
    write_person_table(mort, "participant_mortality__nhanes", evidence_type="person_derived_feature",
                       notes="NCHS public-use Linked Mortality File (NDI linkage), follow-up through 2019-12-31. "
                             "Some decedents' cause/follow-up fields are perturbed in the public-use file; ages < 18 "
                             "are not released.",
                       description="NHANES 2011-2014 linked mortality (public-use LMF 2019)")
    log(f"participants={len(participants)} clinical={len(clinical)} conditions={len(conditions)} "
        f"labs={len(labs)} meds={len(meds)} mortality={len(mort)}")
    return {"participants": participants, "clinical": clinical, "conditions": conditions, "labs": labs,
            "medications": meds, "mortality": mort}



# --------------------------------------------------------------------------------------------
# DATA_AUDIT.md and registry entry (all numbers computed from the tables built above)
# --------------------------------------------------------------------------------------------

LANDING_URL = "https://wwwn.cdc.gov/nchs/nhanes/default.aspx"
DUA_URL = "https://www.cdc.gov/nchs/policy/data-user-agreement.html"
ACCESS_CONDITIONS = (
    "Public-use files, no registration. NCHS Data User Agreement (read 2026-09-23 at " + DUA_URL + "): use the "
    "data for statistical reporting and analysis only; make no attempt to learn the identity of any person; do not "
    "link the data with individually identifiable data from other datasets; do not assess disclosure methods or "
    "attempt re-identification.")


def _pct(n: int, d: int) -> str:
    return f"{n} ({100.0 * n / d:.1f}%)" if d else f"{n} (n/a)"


def _clinical_sources() -> dict[str, str]:
    src = {name: f"{c}.{v}" for c, v, name in YN_ITEMS}
    src.update({name: f"{c}.{v}" for c, v, name, *_ in RANGE_ITEMS})
    src.update({"sbp_mean_mmhg": "BPX.BPXSY1-4", "dbp_mean_mmhg": "BPX.BPXDI1-4", "exam_pulse_60s_bpm": "BPX.BPXPLS",
                "exam_pulse_irregular": "BPX.BPXPULS", "dx_diabetes": "DIQ.DIQ010", "dx_borderline_diabetes": "DIQ.DIQ010",
                "arthritis_type": "MCQ.MCQ195", "smoking_status": "SMQ.SMQ020+SMQ040",
                "phq9_total": "DPQ.DPQ010-090", "phq9_ge10": "DPQ.DPQ010-090", "phq9_n_items": "DPQ.DPQ010-090",
                "phq9_item3_sleep_0to3": "DPQ.DPQ030", "phq9_item4_tired_little_energy_0to3": "DPQ.DPQ040",
                "age_years": "DEMO.RIDAGEYR", "sex": "DEMO.RIAGENDR", "race_ethnicity": "DEMO.RIDRETH3",
                "education_adult20plus": "DEMO.DMDEDUC2", "income_poverty_ratio": "DEMO.INDFMPIR",
                "rx_n_records": "RXQ_RX", "rx_n_distinct_generic": "RXQ_RX.RXDDRUG", "rx_count_reported": "RXQ_RX.RXDCOUNT",
                "mort_eligstat": "LMF.ELIGSTAT", "mort_status": "LMF.MORTSTAT", "mort_ucod_leading": "LMF.UCOD_LEADING",
                "mort_permth_exm": "LMF.PERMTH_EXM", "mort_permth_int": "LMF.PERMTH_INT",
                "pfq061_n_items_answered": "PFQ.PFQ061A-T",
                "healthcare_visits_past_year_cat_0none_1one_2two3_3four9_4ten12_5thirteenplus": "HUQ.HUQ050/HUQ051",
                "any_functional_limitation_pfq": "PFQ.PFQ049+051+054+057+059", "pfq061_n_items_some_difficulty_or_worse": "PFQ.PFQ061A-T"})
    for c in PFQ061_ITEMS:
        src[f"{c.lower()}_difficulty_1none_4unable"] = f"PFQ.{c}"
    for k, name in LAB_WIDE.items():
        comp = next((cp for cp, vs in LAB_SPECS.items() if k in vs), "VITB12")
        src[f"lab_{name}"] = f"{comp}.{k if k != 'VITB12' else 'LBXB12/LBDB12'}"
    return src


SKIP_NOTES = {  # verified against the codebook check items / frequencies
    "DIQ.DIQ070": "asked only if DIQ010 = yes/borderline or DIQ160 = yes (check item DIQ.065)",
    "PFQ.PFQ059": "asked only if PFQ049, PFQ051, PFQ054 and PFQ057 are all 'no' (check item PFQ.058)",
    "HUQ.HUD080": "asked only if hospitalized overnight (HUQ071 = yes; codebook counts match)",
    "MCQ.MCQ195": "asked only if MCQ160A (arthritis) = yes; codebook counts match",
    "LMF.UCOD_LEADING": "defined only for decedents",
    "DPQ.DPQ100": "asked only if >= 1 PHQ-9 item (DPQ010-DPQ090) is 1-3 (see dpq100_skip_check)",
    "HUQ.HUQ050/HUQ051": "HUQ050 (2011-12) and HUQ051 (2013-14) harmonised to the 6 HUQ050 categories",
}


def dpq100_skip_check() -> str:
    """Measured evidence for the DPQ100 skip pattern (both cycles pooled)."""
    d = load_component("DPQ", PHQ9_ITEMS + ["DPQ100"])
    items = d[PHQ9_ITEMS]
    all_zero = items.eq(0).all(axis=1)
    any_pos = items.isin([1.0, 2.0, 3.0]).any(axis=1)
    return (f"asked only if >= 1 PHQ-9 item (DPQ010-DPQ090) is 1-3: missing for "
            f"{int(d.loc[all_zero, 'DPQ100'].isna().sum())} of {int(all_zero.sum())} respondents answering 0 to all "
            f"nine items vs {int(d.loc[any_pos, 'DPQ100'].isna().sum())} of {int(any_pos.sum())} with any item 1-3 "
            "(both cycles)")


def _structural_note(source: str) -> str:
    """Target ages / cycle availability from the codebook (per variable), plus verified skip patterns."""
    if "." not in source:
        return ""
    comp, var = source.split(".", 1)
    notes = []
    if source in SKIP_NOTES:
        notes.append(SKIP_NOTES[source])
    alternatives = [v.split("-")[0].split("+")[0] for v in var.split("/")]
    var = alternatives[0]
    if comp not in COMPONENTS:
        return "; ".join(notes + (["LMF: adults >= 18 only"] if comp == "LMF" else []))
    present = []
    tgt = ""
    for cyc in COMPONENTS[comp]["cycles"]:
        try:
            cb = parse_codebook(comp, cyc)
        except FileNotFoundError:
            continue
        hit = next((v for v in alternatives if v in cb), None)
        if hit:
            present.append(cyc)
            tgt = tgt or cb[hit].get("target", "")
    if present and len(present) < 2:
        notes.append(f"only {CYCLES[present[0]]['label']}")
    m = re.search(r"(\d+) YEARS - (\d+) YEARS", tgt.replace("\n", " "))
    if m and (m.group(1) != "0" or m.group(2) != "150"):
        notes.append(f"target ages {m.group(1)}-{m.group(2)}")
    if comp in LAB_SUBSAMPLE:
        notes.append(LAB_SUBSAMPLE[comp])
    if comp == "PFQ" and var.startswith("PFQ061"):
        notes.append("skip pattern: asked if age >= 60 or 'yes' to PFQ049/PFQ057/PFQ059 (check item PFQ.059A)")
    return "; ".join(notes)


def audit_stats(tables: dict, wear: dict) -> dict:
    p = tables["participants"]
    feats = wear["features"]
    allp = wear["all_participants"]
    days = wear["days"]
    st: dict = {"cycles": {}}
    for cyc in CYCLES:
        pc = p[p["cycle_code"] == cyc]
        ac = allp[allp["cycle_code"] == cyc]
        fc = feats[feats["cycle_code"] == cyc]
        age = pc.set_index("seqn")["age_years"]
        st["cycles"][cyc] = {
            "demo_participants": int(len(pc)),
            "interviewed_only": int((~pc["mec_examined"]).sum()),
            "mec_examined": int(pc["mec_examined"].sum()),
            "paxhd_rows": int(pc["pax_status"].notna().sum()),
            "pam_data_available": int(pc["has_pam_data"].sum()),
            "with_day_rows": int(ac["SEQN"].nunique()),
            "with_minute_rows": wear["qa"].get(cyc, {}).get("participants_with_minute_rows"),
            "passing_rule": int(len(fc)),
            "passing_rule_age18plus": int((fc["seqn"].map(age) >= 18).sum()),
            "passing_rule_age20plus": int((fc["seqn"].map(age) >= 20).sum()),
            "valid_days_ge": {k: int((ac["n_valid_days"] >= k).sum()) for k in range(1, 8)},
            "person_days": int(len(days[days["cycle_code"] == cyc])),
            "valid_person_days": int(days[(days["cycle_code"] == cyc) & days["valid_day"]].shape[0]),
        }
    return st


def write_audit(tables: dict, wear: dict) -> Path:
    from ..wearables.nhanes_features import MIMS_ACTIVE_THRESHOLD, MIMS_THRESHOLD_SOURCE
    MIMS_THRESHOLD_TEXT = f"{MIMS_ACTIVE_THRESHOLD} MIMS/min; {MIMS_THRESHOLD_SOURCE}"
    st = audit_stats(tables, wear)
    p, clin, cond, labs, meds, mort = (tables[k] for k in ["participants", "clinical", "conditions", "labs",
                                                            "medications", "mortality"])
    feats = wear["features"]
    man = load_manifest(SOURCE_ID)["files"]
    qa = wear["qa"]
    age = clin.set_index("participant_id")["age_years"]
    valid_ids = set(feats["participant_id"])
    adult_ids = {i for i in valid_ids if age.get(i, np.nan) >= 18}
    cv = clin[clin["participant_id"].isin(adult_ids)]
    L: list[str] = []
    a = L.append
    tot = {k: sum(c[k] for c in st["cycles"].values() if isinstance(c[k], int)) for k in
           ["demo_participants", "mec_examined", "pam_data_available", "passing_rule", "passing_rule_age18plus"]}
    a("# DATA AUDIT — NHANES 2011-2012 (G) and 2013-2014 (H)\n")
    a(f"_Generated by `measure_it.ingestion.nhanes.write_audit` on {utc_now_iso()}; every count below is computed "
      "from the files retrieved into this directory._\n")
    a("| Field | Value |\n|---|---|")
    rows = [
        ("source_id", SOURCE_ID),
        ("Source (dataset/API name, exact files/endpoints)", "NHANES continuous survey, cycles 2011-2012 (suffix G) and "
         "2013-2014 (suffix H): DEMO, BPX, BMX, MCQ, BPQ, CDQ, DIQ, SLQ, PFQ, DLQ (H only), HSQ, HUQ, DPQ, PAQ, SMQ, HIQ, "
         "INQ, RXQ_RX, CBC, BIOPRO, GHB, GLU, TCHOL, HDL, TRIGLY, VID, VITB12, THYROD (G only), PAXHD, PAXDAY, PAXHR, "
         "PAXMIN; NCHS public-use Linked Mortality Files (2019). See file list below."),
        ("Publishing organization", "CDC National Center for Health Statistics (NCHS)"),
        ("Retrieval date (UTC)", min(v["retrieved_at"] for v in man.values())[:10] + " to "
         + max(v["retrieved_at"] for v in man.values())[:10]),
        ("Source version / release", "NHANES 2011-2012 and 2013-2014 public-use releases (per-file Last-Modified "
         "headers in MANIFEST.json; PAM files first published Nov 2020, PAXMIN large files last-modified 2022-08-01); "
         "Linked Mortality public-use files with follow-up through 2019-12-31"),
        ("Source update date / cadence", "Closed survey cycles; individual files carry their own Last-Modified dates "
         "(MANIFEST.json; latest " + _latest_last_modified(man) + "). The linked_mortality FTP directory lists the 2019 "
         "linkage as the current public-use LMF (older linkages under archived_files/)."),
        ("License / access conditions", ACCESS_CONDITIONS),
        ("Unit of observation", "survey participant (SEQN); labs/meds/conditions long per participant; accelerometry per "
         "participant-minute (PAXMIN), -hour (PAXHR), -day (PAXDAY)"),
        ("Sample size (actual, as ingested)", f"{tot['demo_participants']} participants ({st['cycles']['G']['demo_participants']} "
         f"G + {st['cycles']['H']['demo_participants']} H); {tot['mec_examined']} MEC-examined; {tot['pam_data_available']} with PAM "
         f"data; {tot['passing_rule']} pass the valid-wear rule ({tot['passing_rule_age18plus']} aged >= 18)"),
        ("Geography (resolution, vintage)", "none released — national probability sample; public files carry no geography "
         "and none is inferred (masked variance PSU/strata are design variables, not places)"),
        ("Person-level?", "yes"), ("Geographic?", "no"), ("Omics?", "no"), ("Wearable?", "yes (wrist accelerometry only; no heart rate)"),
        ("True participant linkage across modalities?", "yes — within NHANES only: questionnaire, exam, lab, accelerometry "
         "and linked-mortality records share SEQN (counts below). No linkage to any other dataset."),
    ]
    for k, v in rows:
        a(f"| {k} | {v} |")
    a("\n## Files / endpoints retrieved\n")
    a("Component XPT files: `https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/{2011|2013}/DataFiles/<FILE>.xpt`; "
      "documentation pages (read for every variable used, saved under `docs/`): the `.htm` next to each file. "
      "PAXMIN_G/H are not served at the DataFiles path (HTTP 404, 1,245-byte HTML body); the NHANES examination data "
      "listing links them to `https://ftp.cdc.gov/pub/NHANES/LargeDataFiles/`, which is where they were downloaded "
      "(segmented HTTP range requests, 16 connections per file, because ftp.cdc.gov served ~0.1 MB/s per connection). "
      "Linked mortality: `https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality/`, layout from "
      "the NCHS `R_ReadInProgramAllSurveys.R` read-in program (also saved).\n")
    a("| file | bytes | sha256 (first 16) | Last-Modified | URL | documentation page |\n|---|---:|---|---|---|---|")
    for fn in sorted(man):
        if fn.startswith("docs/"):
            continue
        e = man[fn]
        doc = man.get(f"docs/{fn[:-4]}.htm", {}).get("url", "") if fn.endswith(".xpt") else ""
        a(f"| {fn} | {e['bytes']:,} | {e['sha256'][:16]} | {e.get('last_modified') or ''} | {e['url']} | {doc} |")
    a(f"\nDocumentation pages saved: {sum(1 for f in man if f.startswith('docs/'))} (`docs/*.htm`, URLs in MANIFEST.json).\n")
    a("Requested but not published for these cycles (probed 2026-09-23, each URL returned HTTP 404 with an HTML body):\n")
    for k, v in NOT_PUBLISHED.items():
        a(f"* `{k}` — {v}")
    a("\nConsequences: no ferritin and no CRP/hs-CRP for 2011-2014 in this ingest; thyroid profile only "
      "for 2011-2012 subsample A; disability items (DLQ) only for 2013-2014; prescription reason-for-use ICD-10-CM codes "
      "(RXDRSC1-3) only in RXQ_RX_H (2013-2014) — RXQ_RX_G has no reason variables.\n")
    a("\n## Key variables\n")
    a("| variable (NHANES) | meaning | units | downstream use |\n|---|---|---|---|")
    kv = [
        ("SEQN", "respondent sequence number", "id", "ONLY join key; participant_id = nhanes:<SEQN>"),
        ("SDMVPSU, SDMVSTRA", "masked variance pseudo-PSU / stratum", "code", "survey variance estimation"),
        ("WTINT2YR, WTMEC2YR", "2-year interview / MEC exam weights", "persons", "weighting; pooled 4-year = WTMEC2YR/2 (wtmec4yr_pooled)"),
        ("RIDAGEYR", "age at screening (80 = 80+ top-coded)", "years", "demographics; adult filter"),
        ("BPXSY1-4, BPXDI1-4", "systolic/diastolic readings", "mm Hg", "mean of available readings; diastolic 0 treated as missing"),
        ("BPXPLS", "60-s radial pulse at MEC exam (30 s x 2)", "beats/min", "exam_pulse_60s_bpm — an exam pulse, NOT wearable heart rate; recorded 0 values set to missing"),
        ("BMXBMI, BMXWAIST", "BMI, waist circumference", "kg/m2, cm", "clinical features"),
        ("DPQ010-DPQ090", "PHQ-9 items (0-3)", "score", "phq9_total (complete items only), phq9_ge10"),
        ("DPQ040", "feeling tired or having little energy (past 2 weeks)", "0-3", "fatigue-symptom proxy (>= 2 listed in conditions table as a proxy)"),
        ("HSD010 / HUQ010", "self-rated general health (MEC CAPI 12+ / household)", "1 excellent - 5 poor", "clinical features"),
        ("SLD010H, SLQ050, SLQ060", "usual weekday sleep hours; told doctor trouble sleeping; told sleep disorder", "h, yes/no", "sleep complaints"),
        ("PFQ049-PFQ090, PFQ061A-T", "functional limitations; difficulty items (1 none - 4 unable; 5 'do not do' -> missing)", "yes/no, 1-4", "functional status"),
        ("MCQ*, DIQ010/160, BPQ020/080", "'ever told by a doctor/health professional' diagnoses", "yes/no", "participant_conditions (self-report)"),
        ("RXDDRUG, RXDRSC1-3", "generic drug name; ICD-10-CM reason-for-use codes (2013-14)", "text, code", "participant_medications; conditions labelled 'condition inferred from prescription reason-for-use code'"),
        ("PAXPREDM", "NHANES wake/sleep/non-wear classifier per minute (1/2/3/4 unknown)", "code", "wear masks; algorithm-estimated sleep proxy"),
        ("PAXMTSM", "MIMS triaxial per minute (-0.01 = not computable)", "MIMS", "all minute-level activity features"),
        ("PAXQFM", "count of QC flags in the minute (>0 = invalid)", "count", "minute validity"),
        ("PAXDAY: PAXTMD, PAXWWMD, PAXSWMD, PAXMTSD", "minutes with data; valid wake-/sleep-wear minutes; MIMS day sum over valid minutes", "min, MIMS", "valid-day rule; daily MIMS features"),
        ("PAXHD: PAXSTS, PAXFTIME", "PAM data status; first-data timestamp", "code, HH:MM:SS", "PAM availability; clock alignment of minutes"),
        ("LMF: ELIGSTAT, MORTSTAT, PERMTH_EXM", "linkage eligibility, vital status, months from exam to death/censoring", "code, months", "optional mortality outcome (follow-up to 2019-12-31)"),
    ]
    for r_ in kv:
        a("| " + " | ".join(r_) + " |")
    a("\nFull variable-level documentation (labels, targets, value codes) is in the saved `docs/*.htm` pages.\n")
    from ..wearables.nhanes_features import FEATURE_DEFINITIONS
    a("\n### Derived wearable features (participant_wearable_features__nhanes)\n")
    a(f"Valid-day rule: {wear_rule_text()}. Active/inactive cut-point: {MIMS_THRESHOLD_TEXT}. "
      "Wear minute = PAXPREDM in {1, 2}, PAXQFM == 0, PAXMTSM >= 0. Minute MIMS are stored as float32 in the interim "
      "parquet and rounded back to their published 3-decimal values before thresholding, so a minute at exactly "
      f"{MIMS_ACTIVE_THRESHOLD} counts as active. No heart-rate or HRV features exist or are derived.\n")
    a("| feature | definition | units |\n|---|---|---|")
    for k, (d_, u_) in FEATURE_DEFINITIONS.items():
        a(f"| {k} | {d_} | {u_} |")
    a("\n## Sample sizes (measured)\n")
    a("| count | 2011-2012 (G) | 2013-2014 (H) | total |\n|---|---:|---:|---:|")
    labels = [("demo_participants", "participants in DEMO"), ("interviewed_only", "interviewed only (no MEC exam)"),
              ("mec_examined", "MEC examined"), ("paxhd_rows", "in PAM header file (PAXHD; eligible ages 6+)"),
              ("pam_data_available", "PAM data available (PAXSTS = 1)"), ("with_day_rows", "with any PAXDAY rows"),
              ("with_minute_rows", "with any PAXMIN rows"), ("person_days", "person-days in PAXDAY"),
              ("valid_person_days", "valid person-days (rule below)"), ("passing_rule", "pass valid-wear rule (>= 4 valid days)"),
              ("passing_rule_age18plus", "pass rule and aged >= 18"), ("passing_rule_age20plus", "pass rule and aged >= 20")]
    for key, lab in labels:
        g, h = st["cycles"]["G"][key], st["cycles"]["H"][key]
        t_ = (g + h) if isinstance(g, int) and isinstance(h, int) else "n/a"
        a(f"| {lab} | {g if g is not None else 'n/a'} | {h if h is not None else 'n/a'} | {t_} |")
    a(f"\nValid-day rule: {wear_rule_text()}. Sensitivity (participants with PAM day data having >= k valid days):\n")
    a("| k valid days | " + " | ".join(str(k) for k in range(1, 8)) + " |\n|---|" + "---:|" * 7)
    for cyc in CYCLES:
        a(f"| {CYCLES[cyc]['label']} | " + " | ".join(str(st['cycles'][cyc]['valid_days_ge'][k]) for k in range(1, 8)) + " |")
    a("\nLinkage counts on SEQN (all participants): "
      f"clinical rows {len(clin)}; with >= 1 lab value {labs['participant_id'].nunique()}; with >= 1 prescription record "
      f"{meds['participant_id'].nunique()}; with >= 1 condition row {cond['participant_id'].nunique()}; in LMF "
      f"{len(mort)} (eligible {int((mort['mort_eligstat'] == 'eligible').sum())}, deceased {int((mort['mort_status'] == 1).sum())}). "
      f"Among the {len(adult_ids)} wearable-valid adults: CBC hemoglobin available for "
      f"{int(cv['lab_hemoglobin_g_dl'].notna().sum())}; PHQ-9 total available for "
      f"{int(cv['phq9_total'].notna().sum())}; LMF-eligible {int((cv['mort_eligstat'] == 'eligible').sum())} with "
      f"{int((cv['mort_status'] == 1).sum())} deaths by 2019-12-31.\n")
    for cyc in CYCLES:
        cnt = parse_codebook("PAXMIN", cyc).get("PAXDAYM", {}).get("counts", {})
        expected = sum(v for k, v in cnt.items() if k.isdigit())
        f = xpt_path("PAXMIN", cyc)
        got = xport_layout(f)["n_rows"] if f.exists() else None
        cbm = parse_codebook("PAXMIN", cyc)
        sp_ = interim_stats(cyc)
        pred_cb = {k: v for k, v in cbm.get("PAXPREDM", {}).get("counts", {}).items() if k.isdigit()}
        pred_ok = bool(sp_) and all(sp_["pred_counts"].get(k) == v for k, v in pred_cb.items())
        mims_cb = cbm.get("PAXMTSM", {}).get("counts", {}).get("-0.01")
        mims_ok = bool(sp_) and sp_.get("n_mims_not_computable") == mims_cb
        a(f"* PAXMIN_{cyc} integrity: {got if got is not None else 'file not downloaded'} rows in the downloaded file vs "
          f"{expected:,} minute records in the codebook frequency table (PAXDAYM) — "
          f"{'match' if got == expected else 'MISMATCH'}; PAXPREDM category counts vs codebook — "
          f"{'all match' if pred_ok else 'MISMATCH or not converted'}; PAXMTSM = -0.01 count {sp_.get('n_mims_not_computable')} "
          f"vs codebook {mims_cb} — {'match' if mims_ok else 'MISMATCH or not converted'}. File size equals the server "
          f"Content-Length ({man.get(f'PAXMIN_{cyc}.xpt', {}).get('bytes', 'n/a')} bytes); sha256 in MANIFEST.json.")
    for cyc in CYCLES:
        q = qa.get(cyc, {})
        sp = interim_stats(cyc)
        if sp:
            a(f"* PAXMIN_{cyc}: {sp['n_rows']:,} minute rows ({sp['n_participants']} participants); MIMS not computable "
              f"(-0.01) in {sp['n_mims_not_computable']:,} minutes; QC-flagged minutes {sp['n_qc_flagged']:,}; classifier "
              f"counts {sp['pred_counts']} (1 wake, 2 sleep, 3 non-wear, 4 unknown). Day-start alignment check: "
              f"{q.get('day_starts_checked', 0) - q.get('day_starts_misaligned', 0)}/{q.get('day_starts_checked', 0)} day "
              f"starts at the expected minute. Minute-derived day sums vs PAXDAY PAXMTSD on valid days: "
              f"{q.get('day_sum_check', 'n/a')}.")
    a("\n## Missingness\n")
    a(f"Measured in the wearable-valid adult set (n = {len(cv)}: participants passing the valid-wear rule and aged "
      ">= 18). Structural missingness (age targets, cycle-only components, subsamples, skip patterns) is noted from "
      "the codebooks; 'refused' / 'don't know' are counted as missing.\n")
    a("| variable | source | n non-missing | missing (n, %) | structural reason |\n|---|---|---:|---:|---|")
    SKIP_NOTES["DPQ.DPQ100"] = dpq100_skip_check()
    srcmap = _clinical_sources()
    skip = {"participant_id", "seqn", "cycle", "cycle_code", "sdmvpsu", "sdmvstra", "wtint2yr", "wtmec2yr",
            "has_pam_data", "mec_examined", "age_topcoded_80", "record_key"} | set(c for c in clin.columns if c in (
            "data_layer", "source_name", "source_record_id", "source_version", "retrieved_at",
            "source_geographic_resolution", "evidence_type", "evidence_level", "provenance_notes"))
    for c in clin.columns:
        if c in skip:
            continue
        n_ok = int(cv[c].notna().sum())
        src = srcmap.get(c, "")
        a(f"| {c} | {src} | {n_ok} | {_pct(len(cv) - n_ok, len(cv))} | {_structural_note(src)} |")
    a("\nWearable features (participant_wearable_features__nhanes, all passing participants, n = "
      f"{len(feats)}):\n")
    a("| feature | n non-missing | missing (n, %) |\n|---|---:|---:|")
    for c in feats.columns:
        if feats[c].dtype.kind in "fi" and c not in ("seqn",):
            n_ok = int(feats[c].notna().sum())
            a(f"| {c} | {n_ok} | {_pct(len(feats) - n_ok, len(feats))} |")
    a("\n### Descriptive sanity checks (unweighted, wearable-valid adults; not inferential)\n")
    from scipy.stats import spearmanr
    fa = feats.merge(clin[["participant_id", "age_years", "phq9_item4_tired_little_energy_0to3",
                           "self_rated_health_mec_1excellent_5poor"]], on="participant_id", how="inner")
    fa = fa[fa["age_years"] >= 18]
    a("| feature | Spearman rho with age | n |\n|---|---:|---:|")
    for col in ["mean_daily_mims", "sedentary_fraction", "relative_amplitude", "interdaily_stability",
                "intradaily_variability", "astp", "mean_sleep_wear_min"]:
        if col in fa.columns and fa[col].notna().sum() > 10:
            sub = fa[["age_years", col]].dropna()
            rho = spearmanr(sub["age_years"], sub[col]).statistic
            a(f"| {col} | {rho:.3f} | {len(sub)} |")
    a("\nMedian mean_daily_mims by PHQ-9 item 4 (tired / little energy, past 2 weeks):\n")
    a("| DPQ040 | n | median mean_daily_mims | median sedentary_fraction |\n|---|---:|---:|---:|")
    lab4 = {0.0: "0 not at all", 1.0: "1 several days", 2.0: "2 more than half the days", 3.0: "3 nearly every day"}
    for k, lab_ in lab4.items():
        sub = fa[fa["phq9_item4_tired_little_energy_0to3"] == k]
        sf = f"{sub['sedentary_fraction'].median():.3f}" if "sedentary_fraction" in sub.columns and len(sub) else "n/a"
        a(f"| {lab_} | {len(sub)} | {sub['mean_daily_mims'].median():.0f} | {sf} |")
    a("\n### Conditions among wearable-valid adults\n")
    cva = cond[cond["participant_id"].isin(adult_ids)]
    sr = cva[~cva["evidence_basis"].eq(RX_BASIS)].groupby(["condition_label", "evidence_basis"])["participant_id"].nunique()
    a("| condition | evidence basis | participants |\n|---|---|---:|")
    for (lab_, basis), n in sr.sort_values(ascending=False).items():
        a(f"| {lab_} | {basis} | {n} |")
    rxg = cva[cva["evidence_basis"].eq(RX_BASIS) & cva["rx_target_group"].ne("")].groupby("rx_target_group")["participant_id"].nunique()
    h_adults = int(cv["cycle_code"].eq("H").sum())
    a(f"\nTarget groups from prescription reason-for-use ICD-10-CM codes ({RX_BASIS}; 2013-2014 only, "
      f"denominator {h_adults} wearable-valid adults in 2013-2014):\n")
    a("| group | ICD-10-CM prefixes | participants |\n|---|---|---:|")
    for grp, pref in RX_TARGET_GROUPS.items():
        a(f"| {grp} | {', '.join(pref)} | {int(rxg.get(grp, 0))} |")
    allrx = cond[cond["evidence_basis"].eq(RX_BASIS)]
    absent = [grp for grp, pref in RX_TARGET_GROUPS.items()
              if not allrx["icd10cm_code"].astype(str).str.startswith(tuple(pref)).any()]
    a(f"\nGroups with no matching reason code among all {allrx['participant_id'].nunique()} 2013-2014 participants "
      f"with a coded reason: {', '.join(absent) if absent else 'none'}. These counts identify people taking a "
      "prescription for that reason; they are not prevalence estimates.\n")
    a("\n## Linkage strategy\n")
    rng = p.groupby("cycle_code")["seqn"].agg(["min", "max"])
    overlap = len(set(p.loc[p["cycle_code"] == "G", "seqn"]) & set(p.loc[p["cycle_code"] == "H", "seqn"]))
    a(f"* All NHANES tables are joined on SEQN only; `participant_id = nhanes:<SEQN>`. Measured SEQN ranges: G "
      f"{rng.loc['G', 'min']}-{rng.loc['G', 'max']}, H {rng.loc['H', 'min']}-{rng.loc['H', 'max']}; SEQNs shared by the "
      f"two cycles: {overlap}, so the pooled id is unique.\n"
      "* Questionnaire, examination, laboratory, accelerometry and the NCHS linked-mortality file are true within-person "
      "linkages (same SEQN). Counts per modality are in the tables above.\n"
      "* NOT joinable: nothing in NHANES public files links to any geographic, facility, omics or other person-level "
      "dataset in this project; public files contain no geography below the nation. Any geographic use downstream is "
      "ecological and must never treat NHANES participants as residents of a place.\n"
      "* Survey design: SDMVPSU, SDMVSTRA and WTMEC2YR are kept; pooling G and H requires WTMEC2YR/2 "
      "(`wtmec4yr_pooled`). Fasting-subsample labs (GLU, TRIGLY) need WTSAF2YR and THYROD needs WTSA2YR (not carried in "
      "the wide table). The PAM analyses in NHANES use the MEC exam weight; the valid-wear subset is not separately "
      "re-weighted here, so weighted estimates on that subset carry non-wear selection.\n")
    a("\n## Limitations and caveats\n")
    a("* No heart-rate, HRV, ECG or PPG stream exists in NHANES 2011-2014; `exam_pulse_60s_bpm` (BPXPLS) is a single "
      "60-second pulse count taken at the MEC exam.\n"
      f"* Activity intensity uses a single MIMS cut-point (>= {MIMS_THRESHOLD_TEXT}) — an active/inactive split, not "
      "posture-based sedentary behaviour and not MVPA.\n"
      "* Sleep features are an algorithm-estimated sleep proxy from NHANES's own wake/sleep/non-wear classifier "
      "(PAXPREDM), not validated sleep measurement; the device's idle-sleep mode reduces low-motion information.\n"
      "* Diagnoses are self-reported ('ever told by a doctor'); rx reason codes are participant-reported reasons coded "
      "by NCHS (2013-2014 only) and identify treated conditions, not prevalence. Invisible-illness conditions central "
      "to this project (ME/CFS, Long COVID, POTS, EDS) are not asked about and are essentially absent from reason codes; "
      "fatigue is available only as the PHQ-9 item DPQ040 (a symptom proxy).\n"
      "* Refused/don't-know answers are treated as missing. Age is top-coded at 80. Diastolic readings of 0 are treated "
      "as missing. Below-LOD lab values carry the NCHS fill value (LLOD/sqrt(2)); 2013-2014 vitamin B12 is Deming-adjusted "
      "(LBDB12) and 2011-2012 is unadjusted (LBXB12).\n"
      "* Participants aged 6-17 are included in the tables (PAM eligible from age 6); adult-only analyses must filter.\n"
      "* The public-use linked mortality file perturbs some follow-up/cause fields for confidentiality and excludes "
      "participants under 18; follow-up ends 2019-12-31.\n")
    a("\n## Processed outputs\n")
    from ..store import processed_path
    for name in ["participants__nhanes", "participant_clinical_features__nhanes", "participant_conditions__nhanes",
                 "participant_labs__nhanes", "participant_medications__nhanes", "participant_mortality__nhanes",
                 "participant_wearable_features__nhanes", "participant_wearable_daily__nhanes"]:
        meta = processed_path(name).with_suffix(".meta.json")
        rows_ = json.loads(meta.read_text())["rows"] if meta.exists() else "not written"
        a(f"* `data/processed/{name}.parquet` — {rows_} rows")
    a("\nIntermediate: `data/interim/nhanes_2011_2014/paxmin_{G,H}.parquet` (compact minute tables) with `.stats.json`; "
      "`wearable_qa.json` in this directory.\n")
    a("\n## Reproduce\n")
    a("`uv run python -m measure_it.ingestion.nhanes` (downloads are cached; `--skip-large` omits PAXMIN), or the "
      "wearable step alone: `uv run python -m measure_it.wearables.nhanes_features`. Tests: "
      "`uv run pytest tests/test_nhanes.py tests/test_circadian.py`.\n")
    path = raw_dir(SOURCE_ID) / "DATA_AUDIT.md"
    # keep the extended-laboratory section written by measure_it.ingestion.nhanes_labs_extended (if built)
    ext_section = raw_dir(SOURCE_ID) / "labs_extended_audit_section.md"
    if ext_section.exists():
        L.append("\n" + ext_section.read_text())
    path.write_text("\n".join(L) + "\n")
    return path


def wear_rule_text() -> str:
    from ..wearables.nhanes_features import RULE_TEXT
    return RULE_TEXT


def interim_stats(cycle: str) -> dict:
    from ..wearables.nhanes_features import paxmin_parquet_path
    sp = paxmin_parquet_path(cycle).with_suffix(".stats.json")
    return json.loads(sp.read_text()) if sp.exists() else {}


def _latest_last_modified(man: dict) -> str:
    from email.utils import parsedate_to_datetime
    dates = [parsedate_to_datetime(v["last_modified"]) for k, v in man.items()
             if k.endswith((".xpt", ".dat")) and v.get("last_modified")]
    return max(dates).date().isoformat() + " (latest Last-Modified among retrieved data files)" if dates else UNKNOWN


def write_registry(tables: dict, wear: dict) -> Path:
    st = audit_stats(tables, wear)
    man = load_manifest(SOURCE_ID)["files"]
    minute_ok = all(isinstance(wear["qa"].get(c, {}).get("participants_featurised"), int) for c in CYCLES)
    entry = {
        "source_id": SOURCE_ID,
        "name": "NHANES 2011-2012 (G) and 2013-2014 (H)",
        "publisher": "CDC National Center for Health Statistics (NCHS)",
        "landing_url": LANDING_URL,
        "access_urls": ["https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2011/DataFiles/",
                        "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2013/DataFiles/",
                        "https://ftp.cdc.gov/pub/NHANES/LargeDataFiles/",
                        "https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality/"],
        "license": "U.S. government public-use data; NCHS Data User Agreement applies",
        "access_conditions": ACCESS_CONDITIONS,
        "retrieved_at": max(v["retrieved_at"] for v in man.values()),
        "source_version": "NHANES 2011-2012 and 2013-2014 public releases (PAM released Nov 2020; PAXMIN large files "
                          "last-modified 2022-08-01); public-use LMF 2019",
        "update_date": _latest_last_modified(man),
        "data_layer": "person",
        "unit_of_observation": "participant (SEQN); participant-day and participant-minute for accelerometry",
        "sample_size": {
            "participants": sum(c["demo_participants"] for c in st["cycles"].values()),
            "per_cycle": {CYCLES[k]["label"]: v["demo_participants"] for k, v in st["cycles"].items()},
            "mec_examined": sum(c["mec_examined"] for c in st["cycles"].values()),
            "pam_data_available": sum(c["pam_data_available"] for c in st["cycles"].values()),
            "passing_valid_wear_rule": sum(c["passing_rule"] for c in st["cycles"].values()),
            "passing_valid_wear_rule_age18plus": sum(c["passing_rule_age18plus"] for c in st["cycles"].values()),
        },
        "geographic_resolution": "none",
        "person_level": True, "geographic": False, "omics": False, "wearable": True,
        "participant_linkage": "SEQN joins all NHANES components and the NCHS linked mortality file within NHANES",
        "true_participant_linkage_across_modalities": True,
        "status": "ingested" if minute_ok else "partial",
        "processed_outputs": ["participants__nhanes", "participant_clinical_features__nhanes",
                              "participant_conditions__nhanes", "participant_labs__nhanes",
                              "participant_medications__nhanes", "participant_mortality__nhanes",
                              "participant_wearable_features__nhanes", "participant_wearable_daily__nhanes"],
        "audit": f"data/raw/{SOURCE_ID}/DATA_AUDIT.md",
        "ingestion_module": "measure_it.ingestion.nhanes (+ measure_it.wearables.nhanes_features, measure_it.wearables.circadian)",
        "limitations": [
            "No heart-rate/HRV stream; BPXPLS is a single 60-s exam pulse",
            "Diagnoses are self-reported; rx reason-for-use ICD-10-CM codes exist only for 2013-2014",
            "ME/CFS, Long COVID, POTS and EDS are not measured; fatigue only via PHQ-9 item DPQ040 (proxy)",
            "No ferritin or CRP/hs-CRP in 2011-2014; thyroid profile 2011-2012 subsample only",
            "Single MIMS active/inactive cut-point (10.558); sleep features are classifier-based proxies",
            "No geography in public files; any geographic use is ecological",
        ],
    }
    if not minute_ok:
        entry["status_reason"] = "PAXMIN minute files not processed; minute-derived wearable features missing"
    from .nhanes_labs_extended import registry_extension  # no-op until the extended lab tables exist
    return write_registry_entry(registry_extension(entry))


def _cli() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--download-large-only", action="store_true",
                    help="only fetch PAXMIN_G/H (segmented range download)")
    ap.add_argument("--skip-large", action="store_true", help="do not fetch or process PAXMIN files")
    ap.add_argument("--connections", type=int, default=16)
    args = ap.parse_args()
    if args.download_large_only:
        for cyc in "GH":
            download_large_ranged(xpt_url("PAXMIN", cyc), f"PAXMIN_{cyc}.xpt", n_connections=args.connections,
                                  log=lambda m: print(m, flush=True))
        return
    run(include_large=not args.skip_large)


def run(include_large: bool = True) -> None:
    """Fetch (cached), build clinical tables, build wearable features, write audit + registry entry."""
    fetch_all(include_large=include_large)
    tables = build_clinical_tables()
    from ..wearables import nhanes_features
    wear = nhanes_features.build(use_minute=include_large)
    write_audit(tables, wear)
    write_registry(tables, wear)


if __name__ == "__main__":
    _cli()
