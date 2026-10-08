"""NPPES NPI downloadable file (V.2) + NUCC provider taxonomy -> provider/facility layer.

source_ids: ``nppes`` and ``nucc_taxonomy`` (SPEC section F).

Pipeline
--------
1. Read the CMS NPI_Files.html listing and pick the current *full replacement
   monthly* file. Since 2026-03-03 CMS publishes only Version 2 (V.2).
2. Download the monthly ZIP (~1.1 GB) and the monthly deactivation report via
   ``download_file``; extract the ~11.7 GB ``npidata_pfile`` CSV and the
   secondary-practice-location file next to it.
3. Download the current NUCC taxonomy CSV and write ``nucc_taxonomy``.
4. Stream the NPPES CSV with DuckDB (selected columns only), drop deactivated
   NPIs and non-US practice locations (counts kept for the audit), unpivot the 15
   taxonomy slots, and match codes to ``configs/specialty_groups.yaml``.
5. Geocode the practice ZIP5 to the ZCTA internal point + 2024 county with
   ``crosswalk.zip_to_geo`` and write:
     providers                   one row per NPI in any specialty/facility group
     provider_specialty_groups   one row per (NPI, group) with primary flag
     provider_density_county     2024 county x group counts (zeros included)
     provider_taxonomy_counts    national counts per group and per taxonomy code
     nucc_taxonomy               the NUCC code set with group membership
6. Write DATA_AUDIT.md and registry_entry.yaml for both sources from measured stats.

Guardrail: an NPI taxonomy is self-reported and does NOT mean a provider
evaluates or treats Long COVID, ME/CFS or POTS. NPPES addresses are
self-reported and may be billing/administrative addresses. Coordinates are ZIP
centroids (ZCTA internal points), not address geocodes.

Run: ``uv run python -m measure_it.ingestion.nppes``
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import os
import re
import shutil
import subprocess
import tempfile
import zipfile
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd

from .. import http
from ..config import raw_dir, utc_now_iso
from ..download import download_file, load_manifest
from ..facilities import specialties as sp
from ..geography.crosswalk import STATE_FIPS, zip_to_geo
from ..provenance import add_provenance
from ..registry import write_registry_entry
from ..store import read_table, write_table

SOURCE_ID = "nppes"
NUCC_SOURCE_ID = "nucc_taxonomy"
PRODUCER = "ingestion.nppes"

LISTING_URL = "https://download.cms.gov/nppes/NPI_Files.html"
NPPES_BASE = "https://download.cms.gov/nppes/"
NUCC_CSV_PAGE = ("https://www.nucc.org/index.php/code-sets-mainmenu-41/"
                 "provider-taxonomy-mainmenu-40/csv-mainmenu-57")
NUCC_BASE = "https://www.nucc.org"
# Last Version-1 monthly file still hosted by CMS (V.1 was retired 2026-03-03;
# the March 2026 V.1 file returns 404). Used only to document V.1 vs V.2 layout.
V1_REFERENCE_ZIP = "NPPES_Data_Dissemination_February_2026.zip"

CAVEAT = (
    "An NPI taxonomy code is self-reported and does NOT mean the provider evaluates or treats Long COVID, "
    "ME/CFS, POTS or any other invisible illness. NPPES practice addresses are self-reported and may be "
    "billing/administrative addresses. lat/lon are ZIP5-as-ZCTA internal points (ZIP centroids), not address "
    "geocodes. Deactivated NPIs and non-US practice locations are excluded; only the primary practice location "
    "is used."
)

# NPPES column names (identical in V.1 and V.2 headers).
COLS = {
    "npi": "NPI",
    "entity_type": "Entity Type Code",
    "organization_name": "Provider Organization Name (Legal Business Name)",
    "provider_last_name": "Provider Last Name (Legal Name)",
    "provider_first_name": "Provider First Name",
    "credential": "Provider Credential Text",
    "address_line1": "Provider First Line Business Practice Location Address",
    "city": "Provider Business Practice Location Address City Name",
    "state_raw": "Provider Business Practice Location Address State Name",
    "postal_code_raw": "Provider Business Practice Location Address Postal Code",
    "country_code_raw": "Provider Business Practice Location Address Country Code (If outside U.S.)",
    "phone": "Provider Business Practice Location Address Telephone Number",
    "enumeration_date": "Provider Enumeration Date",
    "last_update_date": "Last Update Date",
    "deactivation_date": "NPI Deactivation Date",
    "reactivation_date": "NPI Reactivation Date",
    "certification_date": "Certification Date",
    "is_organization_subpart": "Is Organization Subpart",
    "parent_organization_name": "Parent Organization LBN",
}
N_SLOTS = 15
for _i in range(1, N_SLOTS + 1):
    COLS[f"tax_code_{_i}"] = f"Healthcare Provider Taxonomy Code_{_i}"
    COLS[f"tax_switch_{_i}"] = f"Healthcare Provider Primary Taxonomy Switch_{_i}"
DATE_COLS = ["enumeration_date", "last_update_date", "deactivation_date", "reactivation_date", "certification_date"]

# ---------------------------------------------------------------- US location rules
US_STATES_DC = {a for f, a in STATE_FIPS.items() if f < "60"}          # 50 states + DC
US_TERRITORIES = {"AS", "GU", "MP", "PR", "VI"}
MILITARY_POSTAL = {"AA", "AE", "AP"}                                    # APO/FPO overseas
FREELY_ASSOCIATED = {"FM", "MH", "PW"}
US_COUNTRY_CODES = {"US", "UM"}   # UM rows are kept only if the state resolves (see classify_us_location)
STATE_NAMES = {
    "ALABAMA": "AL", "ALASKA": "AK", "ARIZONA": "AZ", "ARKANSAS": "AR", "CALIFORNIA": "CA", "COLORADO": "CO",
    "CONNECTICUT": "CT", "DELAWARE": "DE", "DISTRICT OF COLUMBIA": "DC", "FLORIDA": "FL", "GEORGIA": "GA",
    "HAWAII": "HI", "IDAHO": "ID", "ILLINOIS": "IL", "INDIANA": "IN", "IOWA": "IA", "KANSAS": "KS",
    "KENTUCKY": "KY", "LOUISIANA": "LA", "MAINE": "ME", "MARYLAND": "MD", "MASSACHUSETTS": "MA",
    "MICHIGAN": "MI", "MINNESOTA": "MN", "MISSISSIPPI": "MS", "MISSOURI": "MO", "MONTANA": "MT",
    "NEBRASKA": "NE", "NEVADA": "NV", "NEW HAMPSHIRE": "NH", "NEW JERSEY": "NJ", "NEW MEXICO": "NM",
    "NEW YORK": "NY", "NORTH CAROLINA": "NC", "NORTH DAKOTA": "ND", "OHIO": "OH", "OKLAHOMA": "OK",
    "OREGON": "OR", "PENNSYLVANIA": "PA", "RHODE ISLAND": "RI", "SOUTH CAROLINA": "SC", "SOUTH DAKOTA": "SD",
    "TENNESSEE": "TN", "TEXAS": "TX", "UTAH": "UT", "VERMONT": "VT", "VIRGINIA": "VA", "WASHINGTON": "WA",
    "WEST VIRGINIA": "WV", "WISCONSIN": "WI", "WYOMING": "WY", "PUERTO RICO": "PR", "GUAM": "GU",
    "AMERICAN SAMOA": "AS", "NORTHERN MARIANA ISLANDS": "MP", "VIRGIN ISLANDS": "VI",
    "US VIRGIN ISLANDS": "VI", "U S VIRGIN ISLANDS": "VI", "USVI": "VI",
}


def normalize_state(raw) -> str | None:
    """Resolve a self-reported NPPES state field to a USPS code, or None.

    Accepts 2-letter codes, full names (case/punctuation-insensitive) and the
    'CA - CALIFORNIA' / 'CALIFORNIA (CA)' forms seen in the file. Anything else
    (city names, county names, 'N/A', ...) is left unresolved rather than guessed.
    """
    if raw is None or (isinstance(raw, float) and np.isnan(raw)):
        return None
    s = re.sub(r"[.]", "", str(raw).upper()).strip()
    s = re.sub(r"\s+", " ", s)
    valid = US_STATES_DC | US_TERRITORIES | MILITARY_POSTAL | FREELY_ASSOCIATED
    if s in valid:
        return s
    if s in STATE_NAMES:
        return STATE_NAMES[s]
    m = re.match(r"^([A-Z]{2})\s*[-–(]", s)
    if m and m.group(1) in valid:
        return m.group(1)
    m = re.search(r"\(([A-Z]{2})\)", s)
    if m and m.group(1) in valid:
        return m.group(1)
    return None


def classify_us_location(country_code, state_raw) -> tuple[str | None, str]:
    """Return (state_abbr, category) for a practice location.

    Kept categories: 'us_state_dc', 'us_territory'. Dropped: 'dropped_non_us_country',
    'dropped_military_apo_fpo', 'dropped_freely_associated_state', 'dropped_unresolved_state'.
    Country 'UM' (US Minor Outlying Islands) is kept only when the state resolves to a
    US state/DC/territory: the file shows it used for addresses that are otherwise US.
    """
    cc = "" if country_code is None or (isinstance(country_code, float) and np.isnan(country_code)) \
        else str(country_code).strip().upper()
    if cc and cc not in US_COUNTRY_CODES:
        return None, "dropped_non_us_country"
    st = normalize_state(state_raw)
    if st is None:
        return None, "dropped_unresolved_state"
    if st in MILITARY_POSTAL:
        return st, "dropped_military_apo_fpo"
    if st in FREELY_ASSOCIATED:
        return st, "dropped_freely_associated_state"
    if st in US_TERRITORIES:
        return st, "us_territory"
    return st, "us_state_dc"


def zip5(postal) -> str | None:
    """First five digits of a US postal code ('021151234' -> '02115'); None if not 5 leading digits."""
    if postal is None or (isinstance(postal, float) and np.isnan(postal)):
        return None
    m = re.match(r"^\s*(\d{5})", str(postal))
    return m.group(1) if m else None


def display_name(entity_type, first, last, org) -> str | None:
    """Individual: 'FIRST LAST'; organisation: legal business name."""
    def clean(x):
        return "" if x is None or (isinstance(x, float) and np.isnan(x)) else str(x).strip()
    if str(entity_type) == "2":
        return clean(org) or None
    name = " ".join(p for p in (clean(first), clean(last)) if p)
    return name or None


def is_deactivated(deactivation_date, reactivation_date) -> bool:
    """Deactivated unless reactivated on/after the deactivation date (dates as datetime/date or None)."""
    if deactivation_date is None or pd.isna(deactivation_date):
        return False
    if reactivation_date is None or pd.isna(reactivation_date):
        return True
    return reactivation_date < deactivation_date


# ---------------------------------------------------------------- discovery / download
def discover_nppes_release(refresh: bool = False) -> dict:
    """Parse the CMS listing page for the current monthly V.2 file and deactivation report."""
    r = http.get(LISTING_URL, reject_html=False, refresh=refresh)
    r.raise_for_status()
    html = r.text
    monthly_v2 = re.findall(r"href=['\"]\.?/?(NPPES_Data_Dissemination_[A-Za-z]+_\d{4}_V2\.zip)['\"]", html)
    monthly_v1 = re.findall(r"href=['\"]\.?/?(NPPES_Data_Dissemination_[A-Za-z]+_\d{4}\.zip)['\"]", html)
    deact = re.findall(r"href=['\"]\.?/?(NPPES_Deactivated_NPI_Report_\d{6}_V2\.zip)['\"]", html)
    weekly = re.findall(r"href=['\"]\.?/?(NPPES_Data_Dissemination_\d{6}_\d{6}_Weekly_V2\.zip)['\"]", html)
    text = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", html))
    notice = re.search(r"(Effective [^.]*Version 1[^.]*\.[^.]*\.)", text)
    label = re.search(r"(NPPES Data Dissemination V\.2 \([^)]*\) - ZIP format \([^)]*\))", text)
    if not monthly_v2:
        raise RuntimeError(f"no monthly V.2 file linked from {LISTING_URL}; page layout changed?")
    return {
        "listing_url": LISTING_URL,
        "listing_fetched_at": r.fetched_at,
        "monthly_v2_file": monthly_v2[0],
        "monthly_v2_url": NPPES_BASE + monthly_v2[0],
        "monthly_v2_label": label.group(1) if label else None,
        "monthly_v1_files_linked": monthly_v1,
        "deactivation_file": deact[0] if deact else None,
        "deactivation_url": NPPES_BASE + deact[0] if deact else None,
        "weekly_v2_files_linked": weekly,
        "version_notice": notice.group(1) if notice else None,
    }


def discover_nucc_release(refresh: bool = False) -> dict:
    """Parse the NUCC CSV page: the first CSV link is the current version."""
    r = http.get(NUCC_CSV_PAGE, reject_html=False, refresh=refresh)
    r.raise_for_status()
    html = r.text
    links = re.findall(r"href=['\"](/images/stories/CSV/nucc_taxonomy_(\d+)\.csv)['\"]", html)
    if not links:
        raise RuntimeError(f"no NUCC CSV link found on {NUCC_CSV_PAGE}")
    text = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", html))
    ver = re.search(r"current version[^:]*:\s*Version\s+([\d.]+)\s+([\d/]+)", text, re.I)
    path, digits = links[0]
    return {
        "page_url": NUCC_CSV_PAGE,
        "page_fetched_at": r.fetched_at,
        "csv_url": NUCC_BASE + path,
        "csv_file": path.rsplit("/", 1)[-1],
        "version": ver.group(1) if ver else f"{digits[:-1]}.{digits[-1]}",
        "effective": ver.group(2) if ver else None,
        "previous_versions_linked": [p.rsplit("/", 1)[-1] for p, _ in links[1:]],
    }


def extract_members(zip_path: Path, dest: Path) -> dict[str, Path]:
    """Extract the main data CSV, its header, the practice-location CSV and the PDFs (idempotent)."""
    dest.mkdir(parents=True, exist_ok=True)
    out: dict[str, Path] = {}
    with zipfile.ZipFile(zip_path) as zf:
        for info in zf.infolist():
            n = info.filename
            low = n.lower()
            if low.startswith("npidata_pfile") and low.endswith("_fileheader.csv"):
                key = "npidata_header"
            elif low.startswith("npidata_pfile") and low.endswith(".csv"):
                key = "npidata"
            elif low.startswith("pl_pfile") and not low.endswith("_fileheader.csv"):
                key = "pl"
            elif low.endswith(".pdf"):
                key = "readme" if "readme" in low else "codevalues"
            else:
                continue
            target = dest / n
            if not (target.exists() and target.stat().st_size == info.file_size):
                zf.extract(info, dest)
            out[key] = target
    return out


def _nucc_frame() -> tuple[pd.DataFrame, dict]:
    rel = discover_nucc_release()
    path = download_file(rel["csv_url"], NUCC_SOURCE_ID, rel["csv_file"])
    nucc = pd.read_csv(path, dtype=str, encoding="utf-8").fillna("")
    for c in nucc.columns:
        nucc[c] = nucc[c].str.strip()
    rel["rows"] = int(len(nucc))
    rel["path"] = str(path)
    return nucc, rel


# ---------------------------------------------------------------- V.1 vs V.2 layout
class _RangeReader(io.RawIOBase):
    """Seekable file over HTTP byte ranges, so zipfile can read one member of a remote ZIP."""

    def __init__(self, url: str):
        import requests  # file byte-range reads, not API traffic
        http.ensure_online(f"HTTP range reads of {url}")
        self._requests = requests
        self.url, self.pos, self.fetched = url, 0, 0
        h = requests.head(url, timeout=60, headers={"User-Agent": http.USER_AGENT})
        h.raise_for_status()
        self.size = int(h.headers["Content-Length"])

    def seekable(self):
        return True

    def readable(self):
        return True

    def tell(self):
        return self.pos

    def seek(self, off, whence=0):
        self.pos = off if whence == 0 else (self.pos + off if whence == 1 else self.size + off)
        return self.pos

    def read(self, n=-1):
        if n is None or n < 0:
            n = self.size - self.pos
        if n == 0 or self.pos >= self.size:
            return b""
        end = min(self.size, self.pos + n) - 1
        for attempt in range(4):
            try:
                r = self._requests.get(self.url, timeout=120, headers={
                    "Range": f"bytes={self.pos}-{end}", "User-Agent": http.USER_AGENT})
                break
            except self._requests.RequestException:
                if attempt == 3:
                    raise
        if r.status_code != 206:
            raise RuntimeError(f"range request to {self.url} returned HTTP {r.status_code}")
        self.fetched += len(r.content)
        self.pos = end + 1
        return r.content

    def readinto(self, b):
        d = self.read(len(b))
        b[: len(d)] = d
        return len(d)


def _readme_lengths(pdf_bytes: bytes) -> dict[str, int] | None:
    """Parse (section, column) -> max length from a README PDF via pdftotext, if installed."""
    exe = shutil.which("pdftotext")
    if not exe:
        return None
    with tempfile.TemporaryDirectory() as td:
        p = Path(td) / "r.pdf"
        p.write_bytes(pdf_bytes)
        txt = subprocess.run([exe, "-layout", str(p), "-"], capture_output=True, text=True, check=True).stdout
    out, section = {}, None
    for line in txt.splitlines():
        m = re.match(r"^2\.(\d)\s+(.*)", line.strip())
        if m:
            section = m.group(2).strip()
        m = re.match(r"^\s{0,3}([A-Z][A-Za-z0-9 ()\-–./,_]+?)\s{2,}(\d+)(?:\s*\([^)]*\))?\s{2,}(NUMBER|VARCHAR|DATE)",
                     line)
        if m:
            out[f"{section} :: {m.group(1).strip()}"] = int(m.group(2))
            continue
        m = re.match(r"^\s{0,3}([A-Z][A-Za-z0-9 ()\-–./,_]+?)\s{2,}(This is not publicly disseminated)", line)
        if m:
            out[f"{section} :: {m.group(1).strip()}"] = "not publicly disseminated (no length given)"
    return out


def compare_format_versions(members: dict[str, Path], con: duckdb.DuckDBPyConnection, csv_sql: str,
                            refresh: bool = False) -> dict:
    """Document V.1 vs V.2: header columns, README max lengths, and values in the current V.2
    file that exceed the V.1 limits. Best effort; cached in format_version_comparison.json."""
    out_path = raw_dir(SOURCE_ID) / "format_version_comparison.json"
    if out_path.exists() and not refresh:
        return json.loads(out_path.read_text())
    res: dict = {"v1_reference_zip": NPPES_BASE + V1_REFERENCE_ZIP, "computed_at": utc_now_iso()}
    v2_header = members["npidata_header"].read_text().strip().replace('"', "").split(",")
    v2_readme = _readme_lengths(members["readme"].read_bytes()) if "readme" in members else None
    try:
        rr = _RangeReader(NPPES_BASE + V1_REFERENCE_ZIP)
        zf = zipfile.ZipFile(rr)
        ref_dir = raw_dir(SOURCE_ID) / "v1_reference"
        ref_dir.mkdir(exist_ok=True)
        hdr_name = next(n for n in zf.namelist() if n.startswith("npidata_pfile") and n.endswith("_fileheader.csv"))
        rd_name = next(n for n in zf.namelist() if n.lower().endswith(".pdf") and "readme" in n.lower())
        hdr_bytes, rd_bytes = zf.read(hdr_name), zf.read(rd_name)
        (ref_dir / hdr_name).write_bytes(hdr_bytes)
        (ref_dir / rd_name).write_bytes(rd_bytes)
        v1_header = hdr_bytes.decode().strip().replace('"', "").split(",")
        v1_readme = _readme_lengths(rd_bytes)
        res.update({
            "v1_zip_bytes": rr.size, "v1_bytes_fetched_by_range": rr.fetched,
            "v1_members": [i.filename for i in zf.infolist()],
            "v1_header_file": hdr_name, "v1_readme_file": rd_name,
            "v1_header_sha256": hashlib.sha256(hdr_bytes).hexdigest(),
            "v1_readme_sha256": hashlib.sha256(rd_bytes).hexdigest(),
            "v1_n_columns": len(v1_header), "v2_n_columns": len(v2_header),
            "header_columns_identical": v1_header == v2_header,
            "columns_only_in_v1": [c for c in v1_header if c not in v2_header],
            "columns_only_in_v2": [c for c in v2_header if c not in v1_header],
        })
        if v1_readme and v2_readme:
            res["readme_max_length_changes"] = {
                k: {"v1": v1_readme.get(k), "v2": v2_readme.get(k)}
                for k in sorted(set(v1_readme) | set(v2_readme))
                if v1_readme.get(k) != v2_readme.get(k) and "Other Provider Identifier" not in k}
        else:
            res["readme_max_length_changes"] = "pdftotext not installed; README lengths not compared"
    except Exception as exc:  # noqa: BLE001 - documentation step must not break ingestion
        res["v1_error"] = f"{type(exc).__name__}: {exc}"
    # Values in the current V.2 file that would not fit the V.1 field lengths (V.1 README limits).
    lim = {"Provider First Name": 20, "Provider Organization Name (Legal Business Name)": 70,
           "Provider Other First Name": 20, "Provider Other Organization Name": 70,
           "Authorized Official First Name": 20, "Parent Organization LBN": 70}
    sel = ", ".join(f'count(*) FILTER (WHERE length("{c}") > {n}) AS "{c}"' for c, n in lim.items())
    row = con.execute(f"SELECT {sel} FROM {csv_sql}").fetchone()
    res["v2_values_exceeding_v1_limits"] = {c: {"v1_limit": n, "n_values_longer": int(v)}
                                            for (c, n), v in zip(lim.items(), row)}
    res["v2_nonblank_deactivation_reason_code"] = int(con.execute(
        f'SELECT count("NPI Deactivation Reason Code") FROM {csv_sql}').fetchone()[0])
    out_path.write_text(json.dumps(res, indent=2))
    return res


# ---------------------------------------------------------------- core processing
def _q(col: str) -> str:
    return '"' + col.replace('"', '""') + '"'


def _csv_sql(path: Path) -> str:
    p = str(path).replace("'", "''")
    return f"read_csv('{p}', header=true, all_varchar=true, parallel=true)"


def _date_sql(alias: str) -> str:
    return f"try_strptime({alias}, '%m/%d/%Y')::DATE"


def process(members: dict[str, Path], cfg: dict, nucc: pd.DataFrame, threads: int | None = None
            ) -> tuple[dict, dict[str, pd.DataFrame], duckdb.DuckDBPyConnection, str]:
    """Run the DuckDB stage. Returns (stats, frames, connection, csv_sql)."""
    con = duckdb.connect()
    con.execute(f"SET threads={threads or min(16, os.cpu_count() or 4)}")
    csv_sql = _csv_sql(members["npidata"])
    header = members["npidata_header"].read_text().strip().replace('"', "").split(",")
    missing = [c for c in COLS.values() if c not in header]
    if missing:
        raise RuntimeError(f"NPPES header lacks expected columns: {missing}")
    sel = ", ".join((_date_sql(_q(src)) if dst in DATE_COLS else _q(src)) + f" AS {dst}" for dst, src in COLS.items())
    con.execute(f"""CREATE TEMP TABLE raw AS
        SELECT *, (deactivation_date IS NOT NULL
                   AND (reactivation_date IS NULL OR reactivation_date < deactivation_date)) AS deactivated
        FROM (SELECT {sel} FROM {csv_sql})""")

    stats: dict = {}
    stats["rows_total"] = con.execute("SELECT count(*) FROM raw").fetchone()[0]
    stats["npi_unique"] = con.execute("SELECT count(DISTINCT npi) = count(*) FROM raw").fetchone()[0]
    stats["deactivated_excluded"] = con.execute("SELECT count(*) FROM raw WHERE deactivated").fetchone()[0]
    stats["deactivated_with_entity_type_blank"] = con.execute(
        "SELECT count(*) FROM raw WHERE deactivated AND entity_type IS NULL").fetchone()[0]
    stats["reactivated_kept"] = con.execute(
        "SELECT count(*) FROM raw WHERE NOT deactivated AND deactivation_date IS NOT NULL").fetchone()[0]
    stats["reactivation_without_deactivation_date"] = con.execute(
        "SELECT count(*) FROM raw WHERE reactivation_date IS NOT NULL AND deactivation_date IS NULL").fetchone()[0]
    stats["active_rows"] = con.execute("SELECT count(*) FROM raw WHERE NOT deactivated").fetchone()[0]
    stats["active_with_blank_entity_type"] = con.execute(
        "SELECT count(*) FROM raw WHERE NOT deactivated AND entity_type IS NULL").fetchone()[0]
    stats["active_by_entity_type"] = dict(con.execute(
        "SELECT coalesce(entity_type,'blank'), count(*) FROM raw WHERE NOT deactivated GROUP BY 1 ORDER BY 1").fetchall())

    # --- US practice-location filter (classified in Python over distinct (country, state) pairs)
    pairs = con.execute("""SELECT country_code_raw, state_raw, count(*) n FROM raw WHERE NOT deactivated
                           GROUP BY ALL""").df()
    cls = [classify_us_location(c, s) for c, s in zip(pairs["country_code_raw"], pairs["state_raw"])]
    pairs["state"] = [c[0] for c in cls]
    pairs["us_location_category"] = [c[1] for c in cls]
    con.register("loc_pairs", pairs)
    con.execute("""CREATE TEMP TABLE active AS
        SELECT r.* EXCLUDE (deactivated), l.state, l.us_location_category
        FROM raw r JOIN loc_pairs l
          ON r.country_code_raw IS NOT DISTINCT FROM l.country_code_raw AND r.state_raw IS NOT DISTINCT FROM l.state_raw
        WHERE NOT r.deactivated""")
    stats["us_location_categories"] = dict(con.execute(
        "SELECT us_location_category, count(*) FROM active GROUP BY 1 ORDER BY 2 DESC").fetchall())
    stats["dropped_non_us_top_countries"] = dict(con.execute(
        """SELECT country_code_raw, count(*) FROM active WHERE us_location_category='dropped_non_us_country'
           GROUP BY 1 ORDER BY 2 DESC LIMIT 10""").fetchall())
    stats["um_country_code_rows"] = dict(con.execute(
        "SELECT us_location_category, count(*) FROM active WHERE country_code_raw='UM' GROUP BY 1").fetchall())
    con.execute("CREATE TEMP TABLE kept AS SELECT * FROM active WHERE us_location_category IN ('us_state_dc','us_territory')")
    stats["state_resolved_from_nonstandard_text"] = con.execute(
        "SELECT count(*) FROM kept WHERE state <> upper(trim(state_raw))").fetchone()[0]
    stats["kept_us_rows"] = con.execute("SELECT count(*) FROM kept").fetchone()[0]
    stats["kept_by_entity_type"] = dict(con.execute(
        "SELECT entity_type, count(*) FROM kept GROUP BY 1 ORDER BY 1").fetchall())

    # --- taxonomy unpivot (dedupe repeated codes, e.g. same code with licences in several states)
    union = " UNION ALL ".join(
        f"SELECT npi, {i} AS slot, upper(trim(tax_code_{i})) AS code, upper(trim(tax_switch_{i})) AS sw "
        f"FROM kept WHERE tax_code_{i} IS NOT NULL AND trim(tax_code_{i}) <> ''" for i in range(1, N_SLOTS + 1))
    con.execute(f"CREATE TEMP TABLE slots AS {union}")
    stats["taxonomy_slots_filled"] = con.execute("SELECT count(*) FROM slots").fetchone()[0]
    con.execute("""CREATE TEMP TABLE tax AS SELECT npi, code, bool_or(sw='Y') AS y, min(slot) AS first_slot
                   FROM slots GROUP BY npi, code""")
    stats["taxonomy_npi_code_pairs"] = con.execute("SELECT count(*) FROM tax").fetchone()[0]
    stats["taxonomy_duplicate_slots"] = stats["taxonomy_slots_filled"] - stats["taxonomy_npi_code_pairs"]
    stats["distinct_taxonomy_codes"] = con.execute("SELECT count(DISTINCT code) FROM tax").fetchone()[0]
    stats["n_codes_per_npi"] = {int(k): v for k, v in con.execute(
        "SELECT n, count(*) FROM (SELECT npi, count(*) n FROM tax GROUP BY 1) GROUP BY 1 ORDER BY 1").fetchall()}
    stats["kept_without_any_taxonomy"] = con.execute(
        "SELECT count(*) FROM kept k WHERE NOT EXISTS (SELECT 1 FROM tax t WHERE t.npi=k.npi)").fetchone()[0]
    con.execute("""CREATE TEMP TABLE prim AS
        SELECT npi, count(*) AS n_codes, count(*) FILTER (WHERE y) AS n_y,
               arg_min(code, first_slot) FILTER (WHERE y) AS y_code, arg_min(code, first_slot) AS first_code,
               list_sort(list(code)) AS codes
        FROM tax GROUP BY npi""")
    con.execute("""CREATE TEMP TABLE prim2 AS SELECT *,
            CASE WHEN n_y >= 1 THEN y_code ELSE first_code END AS primary_code,
            CASE WHEN n_y = 1 THEN 'primary_switch_Y' WHEN n_y > 1 THEN 'first_Y_of_multiple'
                 WHEN n_codes = 1 THEN 'single_code' ELSE 'first_listed_no_primary_flag' END AS primary_method
        FROM prim; DROP TABLE prim; ALTER TABLE prim2 RENAME TO prim""")
    stats["primary_taxonomy_method"] = dict(con.execute(
        "SELECT primary_method, count(*) FROM prim GROUP BY 1 ORDER BY 2 DESC").fetchall())

    nucc_codes = nucc[["Code"]].rename(columns={"Code": "code"})
    con.register("nucc_codes", nucc_codes)
    stats["codes_not_in_nucc"] = {c: {"npi_code_pairs": int(n), "as_primary": int(p)} for c, n, p in con.execute(
        """SELECT t.code, count(*), count(*) FILTER (WHERE t.code = p.primary_code) FROM tax t JOIN prim p USING (npi)
           WHERE t.code NOT IN (SELECT code FROM nucc_codes) GROUP BY 1 ORDER BY 2 DESC""").fetchall()}

    # --- group membership
    gc = sp.group_code_frame(cfg)
    con.register("gc", gc)
    con.execute("""CREATE TEMP TABLE mem AS
        SELECT t.npi, g.specialty_group, any_value(g.group_type) AS group_type, any_value(g.spec_listed) AS spec_listed,
               bool_or(t.code = p.primary_code) AS is_primary_match,
               array_to_string(list_sort(list_distinct(list(t.code))), '|') AS matched_taxonomy_codes,
               bool_or(g.role = 'physician') AS matched_via_physician_code
        FROM tax t JOIN gc g ON t.code = g.taxonomy_code JOIN prim p ON p.npi = t.npi
        GROUP BY t.npi, g.specialty_group""")
    con.execute("""CREATE TEMP TABLE memagg AS
        SELECT npi,
               array_to_string(list_sort(list(specialty_group)), '|') AS specialty_groups,
               array_to_string(list_sort(list(specialty_group) FILTER (WHERE is_primary_match)), '|') AS primary_specialty_groups,
               array_to_string(list_sort(list(DISTINCT group_type)), '|') AS group_types,
               bool_or(spec_listed) AS in_spec_specialty_group
        FROM mem GROUP BY npi""")

    # --- secondary practice locations (count only; addresses not used)
    pl_sql = _csv_sql(members["pl"])
    con.execute(f"CREATE TEMP TABLE pl AS SELECT NPI AS npi, count(*) AS n_secondary_practice_locations "
                f"FROM {pl_sql} GROUP BY 1")
    stats["pl_rows"] = con.execute(f"SELECT count(*) FROM {pl_sql}").fetchone()[0]
    stats["pl_distinct_npis"] = con.execute("SELECT count(*) FROM pl").fetchone()[0]

    providers = con.execute("""
        SELECT k.npi, k.entity_type, k.provider_last_name, k.provider_first_name, k.organization_name,
               k.credential, k.is_organization_subpart, k.parent_organization_name,
               p.primary_code AS primary_taxonomy_code, p.primary_method AS primary_taxonomy_method,
               p.n_codes AS n_taxonomy_codes, array_to_string(p.codes, '|') AS all_taxonomy_codes,
               m.specialty_groups, m.primary_specialty_groups, m.group_types, m.in_spec_specialty_group,
               k.address_line1, k.city, k.state_raw, k.state, k.us_location_category, k.country_code_raw,
               k.postal_code_raw, k.phone,
               k.enumeration_date, k.last_update_date, k.deactivation_date, k.reactivation_date, k.certification_date,
               coalesce(l.n_secondary_practice_locations, 0) AS n_secondary_practice_locations
        FROM kept k JOIN memagg m USING (npi) JOIN prim p USING (npi) LEFT JOIN pl l USING (npi)
        ORDER BY k.npi""").df()
    membership = con.execute("SELECT * FROM mem ORDER BY npi, specialty_group").df()

    # --- national counts per taxonomy code (included AND considered-but-excluded codes)
    considered = []
    for gid, g in cfg["groups"].items():
        for c in g["codes"]:
            considered.append({"taxonomy_code": c["code"], "config_status": "included", "specialty_group": gid})
        for c in g.get("excluded_considered") or []:
            considered.append({"taxonomy_code": c["code"], "config_status": "excluded_considered",
                               "specialty_group": gid})
    considered = pd.DataFrame(considered)
    con.register("considered_codes", considered[["taxonomy_code"]].drop_duplicates())
    code_counts = con.execute("""
        SELECT c.taxonomy_code,
               count(t.npi) AS n_npis_any,
               count(t.npi) FILTER (WHERE t.code = p.primary_code) AS n_npis_primary,
               count(t.npi) FILTER (WHERE k.entity_type = '1') AS n_individual_any,
               count(t.npi) FILTER (WHERE k.entity_type = '2') AS n_organizations_any
        FROM considered_codes c LEFT JOIN tax t ON t.code = c.taxonomy_code
        LEFT JOIN prim p ON p.npi = t.npi LEFT JOIN kept k ON k.npi = t.npi
        GROUP BY 1 ORDER BY 1""").df()   # ORDER BY: deterministic row order across runs
    joiner = lambda s: "|".join(sorted(set(s)))  # noqa: E731
    status = pd.DataFrame({"considered_in_groups": considered.groupby("taxonomy_code")["specialty_group"].agg(joiner)})
    status["specialty_groups"] = considered[considered["config_status"] == "included"].groupby(
        "taxonomy_code")["specialty_group"].agg(joiner)
    status["config_status"] = np.where(status["specialty_groups"].notna(), "included", "excluded_considered")
    status = status.reset_index()
    code_counts = code_counts.merge(status, on="taxonomy_code", how="left")

    # Field-level facts for the audit.
    stats["kept_last_update_year_quantiles"] = [str(x) for x in con.execute(
        "SELECT quantile_disc(last_update_date, [0.1, 0.25, 0.5, 0.75, 0.9]) FROM kept").fetchone()[0]]
    stats["max_enumeration_date"] = str(con.execute("SELECT max(enumeration_date) FROM raw").fetchone()[0])
    stats["max_last_update_date"] = str(con.execute("SELECT max(last_update_date) FROM raw").fetchone()[0])
    frames = {"providers": providers, "membership": membership, "code_counts": code_counts}
    return stats, frames, con, csv_sql


def deactivation_crosscheck(report_zip: Path, con: duckdb.DuckDBPyConnection) -> dict:
    """Compare the monthly deactivation report with deactivation flags in the full file."""
    with zipfile.ZipFile(report_zip) as zf:
        name = next(n for n in zf.namelist() if n.lower().endswith((".xlsx", ".xls", ".csv")))
        data = zf.read(name)
    if name.lower().endswith(".csv"):
        df = pd.read_csv(io.BytesIO(data), dtype=str, header=None)
    else:
        df = pd.read_excel(io.BytesIO(data), dtype=str, header=None)
    title = str(df.iloc[0, 0])
    df = df[df[0].astype(str).str.fullmatch(r"\d{10}")]
    rep = pd.DataFrame({"npi": df[0].astype(str), "report_date": pd.to_datetime(df[1], format="%m/%d/%Y",
                                                                                 errors="coerce")})
    con.register("deact_report", rep)
    r = con.execute("""
        SELECT count(*) AS n_report,
               count(*) FILTER (WHERE raw.npi IS NULL) AS not_in_full_file,
               count(*) FILTER (WHERE raw.deactivated) AS flagged_deactivated_in_full_file,
               count(*) FILTER (WHERE raw.npi IS NOT NULL AND NOT raw.deactivated) AS active_in_full_file,
               max(d.report_date) AS max_report_date
        FROM deact_report d LEFT JOIN raw ON raw.npi = d.npi""").df().iloc[0]
    return {"report_file": name, "report_title": title, **{k: (str(v) if k == "max_report_date" else int(v))
                                                         for k, v in r.items()}}


# ---------------------------------------------------------------- output tables
def county_fallback() -> pd.Series:
    """ZCTA -> county for ZCTAs whose Census internal point lies outside every 2024 county polygon.

    crosswalk.zip_to_geo leaves these without a county (coastal ZCTAs whose internal point falls in
    water, e.g. 60611 Chicago). Fallback: the county holding the largest land share of the ZCTA in the
    Census 2020 ZCTA-county relationship file, used only if that county exists in the 2024 county set
    (so legacy Connecticut counties are never assigned).
    """
    zc = read_table("zcta_centroids", columns=["zcta", "county_fips"])
    missing = set(zc.loc[zc["county_fips"].isna(), "zcta"])
    xw = read_table("zcta_county_crosswalk", columns=["zcta", "county_fips_2020", "land_share_of_zcta"])
    xw = xw[xw["zcta"].isin(missing)].sort_values("land_share_of_zcta", ascending=False).drop_duplicates("zcta")
    g = read_table("geographies", columns=["geo_level", "county_fips", "ct_legacy"])
    c24 = set(g.loc[(g["geo_level"] == "county") & (~g["ct_legacy"].astype(bool)), "county_fips"])
    xw = xw[xw["county_fips_2020"].isin(c24)]
    return xw.set_index("zcta")["county_fips_2020"]


def membership_resolution(geocode_method: pd.Series) -> pd.Series:
    """source_geographic_resolution for provider_specialty_groups rows: the table's county_fips is derived
    from the ZIP5-as-ZCTA internal point ('zcta_centroid'); rows that were not geocoded carry no location."""
    return pd.Series(np.where(geocode_method == "zip5_as_zcta_internal_point", "zcta_centroid", "none"),
                     index=geocode_method.index)


def connecticut_summary(prov: pd.DataFrame) -> dict:
    """Measured Connecticut vintage facts: CT-address providers, their county vintage, legacy FIPS assigned."""
    ct = prov[prov["state"] == "CT"]
    cf = ct["county_fips"].dropna()
    legacy = {f"090{i:02d}" for i in range(1, 16, 2)}   # 09001..09015 (2020-vintage CT counties)
    return {
        "ct_address_providers": int(len(ct)),
        "ct_geocoded_to_2024_planning_region": int(cf.str.match(r"^091[1-9]0$").sum()),
        "ct_geocoded_to_legacy_county": int(cf.isin(legacy).sum()),
        "ct_geocoded_to_other_state_county": int((~cf.str.startswith("09")).sum()),
        "ct_not_geocoded": int(ct["county_fips"].isna().sum()),
        "any_provider_on_legacy_ct_county": bool(prov["county_fips"].isin(legacy).any()),
    }


def build_outputs(frames: dict, cfg: dict, nucc: pd.DataFrame, meta: dict) -> dict:
    prov = frames["providers"].copy()
    mem = frames["membership"].copy()
    nucc_idx = nucc.set_index("Code")

    prov["entity_type"] = prov["entity_type"].astype(int)
    prov["entity_type_label"] = prov["entity_type"].map({1: "individual", 2: "organization"})
    prov["name"] = [display_name(e, f, ln, o) for e, f, ln, o in zip(
        prov["entity_type"], prov["provider_first_name"], prov["provider_last_name"], prov["organization_name"])]
    prov.loc[prov["entity_type"] == 2, "credential"] = None  # organisations have no provider credential
    for col, src in [("primary_taxonomy_display_name", "Display Name"), ("primary_taxonomy_grouping", "Grouping"),
                     ("primary_taxonomy_classification", "Classification"),
                     ("primary_taxonomy_specialization", "Specialization")]:
        prov[col] = prov["primary_taxonomy_code"].map(nucc_idx[src]).replace("", None)
    prov["primary_specialty_groups"] = prov["primary_specialty_groups"].fillna("")
    prov["zip5"] = prov["postal_code_raw"].map(zip5)
    geo = zip_to_geo(prov["zip5"])
    prov["lat"], prov["lon"] = geo["lat"].values, geo["lon"].values
    prov["county_fips"] = geo["county_fips"].values
    prov["geocode_method"] = geo["geocode_method"].values
    prov["county_assignment_method"] = np.where(prov["county_fips"].notna(), "zcta_point_in_2024_county", None)
    fb = county_fallback()
    need = (prov["geocode_method"] == "zip5_as_zcta_internal_point") & prov["county_fips"].isna()
    prov.loc[need, "county_fips"] = prov.loc[need, "zip5"].map(fb)
    prov.loc[need & prov["county_fips"].notna(), "county_assignment_method"] = "zcta_largest_land_share_rel2020"
    prov["state_fips"] = prov["state"].map({v: k for k, v in STATE_FIPS.items()})
    prov["in_50_states_dc"] = prov["state"].isin(US_STATES_DC)
    geo_state = prov["county_fips"].str[:2]
    prov["geo_state_matches_address_state"] = pd.array(
        np.where(geo_state.notna(), geo_state == prov["state_fips"], None), dtype="boolean")
    prov["not_geocoded_reason"] = np.select(
        [prov["geocode_method"] != "not_geocoded", prov["zip5"].isna()],
        ["", "no_5_digit_zip"], default="zip5_not_a_2024_zcta")
    prov["phone"] = prov["phone"].str.replace(r"\D", "", regex=True).replace("", None)
    for c in ["provider_last_name", "provider_first_name", "organization_name", "credential", "address_line1",
              "city", "parent_organization_name"]:
        prov[c] = prov[c].where(prov[c].isna(), prov[c].str.strip()).replace("", None)

    version = meta["source_version"]
    retrieved = meta["retrieved_at"]
    prov_cols = [
        "npi", "entity_type", "entity_type_label", "name", "provider_last_name", "provider_first_name",
        "organization_name", "credential", "is_organization_subpart", "parent_organization_name",
        "primary_taxonomy_code", "primary_taxonomy_display_name", "primary_taxonomy_grouping",
        "primary_taxonomy_classification", "primary_taxonomy_specialization", "primary_taxonomy_method",
        "n_taxonomy_codes", "all_taxonomy_codes", "specialty_groups", "primary_specialty_groups", "group_types",
        "in_spec_specialty_group", "address_line1", "city", "state", "state_raw", "state_fips", "in_50_states_dc",
        "us_location_category", "country_code_raw", "zip5", "phone", "enumeration_date", "last_update_date",
        "deactivation_date", "reactivation_date", "certification_date", "n_secondary_practice_locations",
        "lat", "lon", "county_fips", "geocode_method", "county_assignment_method", "not_geocoded_reason",
        "geo_state_matches_address_state",
    ]
    prov = prov[prov_cols]
    res = np.where(prov["geocode_method"] == "zip5_as_zcta_internal_point", "zcta_centroid", "state")
    providers = add_provenance(
        prov, data_layer="facility", source_name="CMS NPPES NPI Downloadable File (V.2) + NUCC taxonomy",
        source_version=version, retrieved_at=retrieved, evidence_type="provider_registry",
        source_record_id="npi", source_geographic_resolution="zcta_centroid",
        evidence_level="self_reported_registry", provenance_notes=CAVEAT)
    providers["source_geographic_resolution"] = res

    # --- membership (long)
    ent = prov.set_index("npi")[["entity_type", "county_fips", "geocode_method"]]
    mem = mem.join(ent, on="npi")
    membership = add_provenance(
        mem[["npi", "specialty_group", "group_type", "spec_listed", "is_primary_match", "matched_taxonomy_codes",
             "matched_via_physician_code", "entity_type", "county_fips"]],
        data_layer="facility", source_name="CMS NPPES NPI Downloadable File (V.2) + configs/specialty_groups.yaml",
        source_version=version, retrieved_at=retrieved, evidence_type="provider_registry",
        source_record_id=lambda d: d["npi"] + ":" + d["specialty_group"], source_geographic_resolution="none",
        evidence_level="self_reported_registry", provenance_notes=CAVEAT)
    # county_fips here comes from the ZIP5-as-ZCTA internal point, so label it as such (rows without a
    # geocode carry no location in this table).
    membership["source_geographic_resolution"] = membership_resolution(mem["geocode_method"]).values

    # --- county density (all 2024 counties x all groups, zeros kept)
    groups = sp.group_code_frame(cfg)[["specialty_group", "group_label", "group_type", "spec_listed"]].drop_duplicates()
    ind_flag = mem["entity_type"] == 1
    org_flag = mem["entity_type"] == 2
    flags = pd.DataFrame({
        "county_fips": mem["county_fips"], "specialty_group": mem["specialty_group"], "npi": mem["npi"],
        "n_individual_providers": ind_flag,
        "n_individual_providers_primary": ind_flag & mem["is_primary_match"],
        "n_individual_physicians": ind_flag & mem["matched_via_physician_code"],
        "n_organizations": org_flag,
        "n_organizations_primary": org_flag & mem["is_primary_match"],
    })
    count_cols = ["n_individual_providers", "n_individual_providers_primary", "n_individual_physicians",
                  "n_organizations", "n_organizations_primary"]
    agg = flags[flags["county_fips"].notna()].groupby(["county_fips", "specialty_group"])[count_cols].sum()
    agg = agg.astype("int64").reset_index()
    counties = read_table("geographies")
    counties = counties[(counties["geo_level"] == "county") & (~counties["ct_legacy"].astype(bool))][
        ["county_fips", "name", "state_abbr"]].rename(columns={"name": "county_name"})
    grid = counties.merge(groups, how="cross")
    dens = grid.merge(agg, on=["county_fips", "specialty_group"], how="left")
    unmatched_counties = sorted(set(agg["county_fips"]) - set(counties["county_fips"]))
    dens[count_cols] = dens[count_cols].fillna(0).astype("int64")
    dens = dens[["county_fips", "county_name", "state_abbr", "specialty_group", "group_label", "group_type",
                 "spec_listed"] + count_cols]
    density = add_provenance(
        dens, data_layer="facility", source_name="CMS NPPES NPI Downloadable File (V.2), aggregated to 2024 counties",
        source_version=version, retrieved_at=retrieved, evidence_type="provider_registry",
        source_record_id=lambda d: d["county_fips"] + ":" + d["specialty_group"],
        source_geographic_resolution="county", evidence_level="self_reported_registry",
        provenance_notes=("Counts of NPIs whose practice ZIP5 internal point falls in the county; any-taxonomy match "
                          "unless *_primary. Groups overlap: never sum across groups. Denominators (per 100k) are "
                          "added downstream from ACS. " + CAVEAT))

    # --- national counts
    flags["n_npis_any"] = True
    flags["n_npis_primary"] = mem["is_primary_match"]
    flags["n_geocoded"] = mem["county_fips"].notna()
    nat = flags.groupby("specialty_group")[["n_npis_any", "n_npis_primary", "n_individual_providers",
                                            "n_individual_providers_primary", "n_individual_physicians",
                                            "n_organizations", "n_organizations_primary", "n_geocoded"]].sum()
    nat = nat.rename(columns={"n_individual_providers": "n_individual_any",
                              "n_individual_providers_primary": "n_individual_primary",
                              "n_organizations": "n_organizations_any"})
    nat["n_counties_with_any"] = flags.groupby("specialty_group")["county_fips"].nunique()
    grp = groups.merge(nat.reset_index(), on="specialty_group", how="left")
    ncols = [c for c in grp.columns if c.startswith("n_")]
    grp[ncols] = grp[ncols].fillna(0).astype("int64")
    grp["n_not_geocoded"] = grp["n_npis_any"] - grp["n_geocoded"]
    grp["count_level"] = "specialty_group"
    grp["taxonomy_code"] = None
    grp["taxonomy_display_name"] = None
    grp["config_status"] = "group"
    cc = frames["code_counts"].copy()
    cc["count_level"] = "taxonomy_code"
    cc["taxonomy_display_name"] = cc["taxonomy_code"].map(nucc_idx["Display Name"])
    cc = cc.rename(columns={"specialty_groups": "specialty_group"})
    cc["group_label"] = None
    cc["group_type"] = None
    cc["spec_listed"] = None
    counts = pd.concat([grp, cc], ignore_index=True, sort=False)
    int_cols = ["n_npis_any", "n_npis_primary", "n_individual_any", "n_organizations_any", "n_geocoded",
                "n_counties_with_any", "n_individual_primary", "n_individual_physicians", "n_organizations_primary",
                "n_not_geocoded"]
    for c in int_cols:
        counts[c] = pd.to_numeric(counts[c]).astype("Int64")
    counts = counts[["count_level", "specialty_group", "group_label", "group_type", "spec_listed", "taxonomy_code",
                     "taxonomy_display_name", "config_status", "considered_in_groups"] + int_cols]
    counts["spec_listed"] = counts["spec_listed"].astype("boolean")
    counts = add_provenance(
        counts, data_layer="facility", source_name="CMS NPPES NPI Downloadable File (V.2) + NUCC taxonomy",
        source_version=version, retrieved_at=retrieved, evidence_type="provider_registry",
        source_record_id=lambda d: d["count_level"] + ":" + d["taxonomy_code"].fillna(d["specialty_group"]).astype(str),
        source_geographic_resolution="national", evidence_level="self_reported_registry",
        provenance_notes=("US practice locations, active NPIs. taxonomy_code rows cover included codes and codes "
                          "considered but excluded (config_status). " + CAVEAT))
    return {"providers": providers, "membership": membership, "density": density, "counts": counts,
            "unmatched_counties": unmatched_counties, "connecticut": connecticut_summary(prov)}


def build_nucc_table(nucc: pd.DataFrame, cfg: dict, rel: dict, retrieved: str) -> pd.DataFrame:
    c2g = sp.code_to_groups(cfg)
    t = pd.DataFrame({
        "taxonomy_code": nucc["Code"], "grouping": nucc["Grouping"], "classification": nucc["Classification"],
        "specialization": nucc["Specialization"].replace("", None), "display_name": nucc["Display Name"],
        "section": nucc["Section"], "definition": nucc["Definition"].replace("", None),
        "notes": nucc["Notes"].replace("", None),
    })
    t["role"] = [sp.role_from_nucc(g, c, s) for g, c, s in zip(t["grouping"], t["classification"], t["section"])]
    t["specialty_groups"] = t["taxonomy_code"].map(lambda c: "|".join(c2g.get(c, ())))
    t["nucc_version"] = rel["version"]
    return add_provenance(
        t, data_layer="facility", source_name="NUCC Health Care Provider Taxonomy Code Set (CSV)",
        source_version=f"{rel['version']} (effective {rel.get('effective')})", retrieved_at=retrieved,
        evidence_type="ontology_mapping", source_record_id="taxonomy_code", source_geographic_resolution="none",
        provenance_notes="specialty_groups from configs/specialty_groups.yaml (hand-curated)")


# ---------------------------------------------------------------- audit + registry
def _pct(n, d) -> str:
    return f"{n:,} ({(100.0 * n / d if d else 0):.2f}%)"


def _top_ungeocoded_zips(ng: pd.DataFrame, n: int = 10) -> str:
    """'44195: 2,379 (CLEVELAND OH); ...' for the ZIP5s holding the most ungeocoded providers."""
    out = []
    for z, cnt in ng["zip5"].value_counts().head(n).items():
        sub = ng[ng["zip5"] == z]
        city = sub["city"].str.upper().mode()
        out.append(f"{z}: {cnt:,} ({city.iloc[0] if len(city) else '?'} {sub['state'].mode().iloc[0]})")
    return "; ".join(out)


def _missingness(df: pd.DataFrame, cols: list[str], mask=None) -> list[str]:
    d = df if mask is None else df[mask]
    rows = []
    for c in cols:
        s = d[c]
        miss = s.isna() | (s.astype(str).str.strip() == "")
        rows.append(f"| {c} | {len(d):,} | {_pct(int(miss.sum()), len(d))} |")
    return rows


def write_nppes_audit(stats: dict, outs: dict, cfg: dict, release: dict, fmt: dict, deact: dict,
                      manifest: dict, meta: dict) -> Path:
    prov = outs["providers"]
    counts = outs["counts"]
    grp = counts[counts["count_level"] == "specialty_group"]
    codes = counts[counts["count_level"] == "taxonomy_code"]
    n = len(prov)
    ind = prov["entity_type"] == 1
    geocoded = prov["geocode_method"] == "zip5_as_zcta_internal_point"
    terr = ~prov["in_50_states_dc"]
    geo_rate = (f"not geocoded: territories {_pct(int((terr & ~geocoded).sum()), int(terr.sum()))} vs 50 states + DC "
                f"{_pct(int((~terr & ~geocoded).sum()), int((~terr).sum()))}")
    im = prov["all_taxonomy_codes"].str.contains("207R00000X", regex=False)
    im_only = im & (prov["primary_specialty_groups"] == "")
    im_hosp = im & prov["all_taxonomy_codes"].str.contains("208M00000X", regex=False)
    im_other_grp = im & prov["specialty_groups"].str.split("|").map(lambda g: len(set(g) - {"primary_care"}) > 0)
    fb_zips = prov.loc[prov["county_assignment_method"] == "zcta_largest_land_share_rel2020",
                       "zip5"].value_counts().head(8).to_dict()
    im_note = (f"Of {int(im.sum()):,} NPIs carrying 207R00000X Internal Medicine: {_pct(int(im_hosp.sum()), int(im.sum()))} "
               f"also carry 208M00000X Hospitalist, {_pct(int(im_other_grp.sum()), int(im.sum()))} also match another "
               f"group (a subspecialty or facility code), and {_pct(int(im_only.sum()), int(im.sum()))} have a primary "
               "taxonomy outside every group.")
    lines = ["# DATA AUDIT — NPPES NPI Downloadable File (V.2)", "",
             "> **Mandatory caveat.** " + CAVEAT, "",
             "| Field | Value |", "|---|---|",
             f"| source_id | `{SOURCE_ID}` (+ `{NUCC_SOURCE_ID}` for code names; see its own audit) |",
             f"| Source (dataset/API name, exact files/endpoints) | {release['monthly_v2_label'] or release['monthly_v2_file']}: "
             f"`{release['monthly_v2_url']}` → `{meta['csv_name']}`; deactivation report `{release['deactivation_url']}` |",
             "| Publishing organization | Centers for Medicare & Medicaid Services (CMS), NPPES |",
             f"| Retrieval date (UTC) | {meta['retrieved_at']} |",
             f"| Source version / release | {meta['source_version']} (format **Version 2**) |",
             f"| Source update date / cadence | ZIP Last-Modified {manifest.get('last_modified')}; latest enumeration date "
             f"{stats['max_enumeration_date']}, latest last-update date {stats['max_last_update_date']}; full replacement "
             f"monthly, incremental weekly |",
             "| License / access conditions | FOIA-disclosable public data, free download, no registration or DUA. CMS: "
             "issuance of an NPI does not validate licensure or credentials. |",
             "| Unit of observation | NPI (individual provider, entity type 1, or organisation/subpart, entity type 2) |",
             f"| Sample size (actual, as ingested) | {stats['rows_total']:,} NPI rows in file; {stats['active_rows']:,} "
             f"active; {stats['kept_us_rows']:,} active with a US practice location; **{n:,} NPIs in ≥1 specialty/"
             f"facility group** (`providers`) |",
             "| Geography (resolution, vintage) | Self-reported practice address; geocoded to ZIP5-as-ZCTA internal point "
             "(2024 Gazetteer) and the 2024 county containing it |",
             "| Person-level? | no (providers are not patients; no patient data) |",
             "| Geographic? | yes (facility locations; county aggregates) |",
             "| Omics? | no |", "| Wearable? | no |",
             "| True participant linkage across modalities? | not applicable — no participants |", "",
             "## Files / endpoints retrieved", "",
             f"* Listing page `{LISTING_URL}` (fetched {release['listing_fetched_at']}). Page notice: "
             f"\"{release['version_notice']}\"",
             f"* Monthly V.2 files linked: `{release['monthly_v2_file']}`; V.1 monthly files linked: "
             f"{release['monthly_v1_files_linked'] or 'none'}; weekly V.2 files linked: {release['weekly_v2_files_linked']}",
             "* Downloaded (sha256 and byte sizes in `MANIFEST.json`):"]
    for fname, e in load_manifest(SOURCE_ID)["files"].items():
        lines.append(f"  * `{fname}` — {e['bytes']:,} bytes, sha256 `{e['sha256']}`, Last-Modified {e.get('last_modified')}")
    lines += [f"* Extracted from the ZIP into `data/raw/nppes/extracted/`: `{meta['csv_name']}` "
              f"({meta['csv_bytes']:,} bytes), its `_fileheader.csv`, `{meta['pl_name']}` (secondary practice "
              "locations, used only for a per-NPI count), README v.2 and Code Values PDFs.", "",
              "### Format version: is there an NPPES V.2? Yes, and it is now the only format", "",
              "* The CMS listing states that Version 1 of the monthly/weekly files is no longer supported from "
              "03/03/2026 and links only V.2 files. The file used is the **V.2 full replacement monthly file**.",
              f"* V.1 reference for comparison: `{fmt.get('v1_reference_zip')}` (last V.1 monthly still hosted; "
              "`NPPES_Data_Dissemination_March_2026.zip` returned HTTP 404 when probed on 2026-09-23). Only its header "
              f"CSV and README were read, by HTTP range requests ({fmt.get('v1_bytes_fetched_by_range', 'n/a'):,} bytes of "
              f"{fmt.get('v1_zip_bytes', 0):,}); copies in `data/raw/nppes/v1_reference/`."
              if isinstance(fmt.get("v1_bytes_fetched_by_range"), int) else
              f"* V.1 reference could not be read: {fmt.get('v1_error')}"]
    if "header_columns_identical" in fmt:
        lines.append(f"* Data-file header: V.1 {fmt['v1_n_columns']} columns, V.2 {fmt['v2_n_columns']} columns; "
                     f"column names identical and in the same order: **{fmt['header_columns_identical']}**. "
                     "(Both include `Healthcare Provider Taxonomy Group_1..15` and `Certification Date`.)")
    chg = fmt.get("readme_max_length_changes")
    if isinstance(chg, dict):
        lines.append("* README max-length differences (V.1 README → V.2 README; parsed with pdftotext):")
        for k, v in chg.items():
            lines.append(f"  * {k}: {v['v1']} → {v['v2']}")
    elif chg:
        lines.append(f"* README comparison: {chg}")
    if "v2_nonblank_deactivation_reason_code" in fmt:
        lines.append(f"* `NPI Deactivation Reason Code`: non-blank values in the V.2 data file: "
                     f"{fmt['v2_nonblank_deactivation_reason_code']:,}.")
    ex = fmt.get("v2_values_exceeding_v1_limits", {})
    if ex:
        lines.append("* Measured impact in the current V.2 file — values longer than the V.1 limit "
                     "(these would have been truncated in V.1): " +
                     "; ".join(f"{k} > {v['v1_limit']}: {v['n_values_longer']:,}" for k, v in ex.items()) + ".")
    lines += ["", "## Processing and exclusions (all counts measured)", "",
              f"* NPI rows in file: {stats['rows_total']:,} (NPI unique: {stats['npi_unique']}).",
              f"* **Deactivated NPIs excluded**: {stats['deactivated_excluded']:,} (deactivation date set and no "
              f"reactivation on/after it; {stats['deactivated_with_entity_type_blank']:,} of these have all other fields "
              f"blanked by CMS). Reactivated NPIs kept: {stats['reactivated_kept']:,}; rows with a reactivation date but "
              f"no deactivation date (kept): {stats['reactivation_without_deactivation_date']:,}.",
              f"* Cross-check against the monthly deactivation report ({deact.get('report_title')}): "
              f"{deact.get('n_report', 0):,} NPIs listed; {deact.get('flagged_deactivated_in_full_file', 0):,} flagged "
              f"deactivated in the full file, {deact.get('active_in_full_file', 0):,} active in the full file, "
              f"{deact.get('not_in_full_file', 0):,} absent from it (latest deactivation date in the report "
              f"{str(deact.get('max_report_date'))[:10]}; data file `{meta['csv_name']}`, latest last-update date "
              f"{stats['max_last_update_date']}).",
              f"* Active rows by entity type: {stats['active_by_entity_type']} (blank entity type among active: "
              f"{stats['active_with_blank_entity_type']:,}).",
              "* **US practice locations only.** Practice-location category counts among active NPIs:"]
    for k, v in stats["us_location_categories"].items():
        lines.append(f"  * `{k}`: {v:,}")
    lines += [f"  * Top non-US country codes dropped: {stats['dropped_non_us_top_countries']}",
              f"  * Country code `UM` (US Minor Outlying Islands) rows: {stats['um_country_code_rows']} — their state "
              "fields are mostly US states or Puerto Rico; kept only when the state resolves to a US state/DC/territory.",
              f"  * State text resolved from a non-standard form (full name, 'CA - CALIFORNIA', 'P.R.' …): "
              f"{stats['state_resolved_from_nonstandard_text']:,} rows. Unresolvable text is dropped, not guessed.",
              "  * Territories (PR, GU, VI, AS, MP) are kept and flagged `in_50_states_dc = False`.",
              f"* Kept (active, US): {stats['kept_us_rows']:,}; by entity type {stats['kept_by_entity_type']}.",
              f"* Taxonomy slots filled: {stats['taxonomy_slots_filled']:,}; distinct (NPI, code) pairs "
              f"{stats['taxonomy_npi_code_pairs']:,} ({stats['taxonomy_duplicate_slots']:,} repeated slots, e.g. the same "
              f"code listed with licences in several states); distinct codes {stats['distinct_taxonomy_codes']:,}; "
              f"kept NPIs with no taxonomy: {stats['kept_without_any_taxonomy']:,}.",
              f"* Taxonomy codes per NPI: {stats['n_codes_per_npi']}.",
              f"* Primary taxonomy selection: {stats['primary_taxonomy_method']} (every kept NPI has exactly one slot "
              "with Primary Taxonomy Switch = Y when the method is `primary_switch_Y`).",
              f"* Codes in NPPES not present in NUCC {meta['nucc_version']}: {stats['codes_not_in_nucc']} (display "
              "names left null for these).",
              f"* Secondary practice locations file: {stats['pl_rows']:,} rows for {stats['pl_distinct_npis']:,} NPIs; "
              "used only for `n_secondary_practice_locations` — providers are placed at their PRIMARY practice address.",
              ""]
    lines += ["## Specialty and facility groups", "",
              "Defined in `configs/specialty_groups.yaml` (hand-chosen NUCC codes with verbatim NUCC display names; "
              "checked against the NUCC file by `tests/test_nppes.py`). A provider matches a group if ANY of its up to 15 "
              "taxonomy codes is listed; `primary` = its primary taxonomy code is listed. Groups overlap (never sum "
              "across groups).", "",
              "| group | type | SPEC | codes | NPIs any | NPIs primary | individuals any | individual physicians | "
              "organisations any | not geocoded | counties with ≥1 |", "|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in grp.itertuples():
        ncodes = len(cfg["groups"][r.specialty_group]["codes"])
        lines.append(f"| {r.specialty_group} | {r.group_type} | {r.spec_listed} | {ncodes} | {r.n_npis_any:,} | "
                     f"{r.n_npis_primary:,} | {r.n_individual_any:,} | {r.n_individual_physicians:,} | "
                     f"{r.n_organizations_any:,} | {r.n_not_geocoded:,} | {r.n_counties_with_any:,} |")
    lines += ["", "### Ambiguous choices (codes considered and excluded, with measured national NPI counts)", "",
              "| group | code | NUCC display name | NPIs any | as primary | reason excluded |", "|---|---|---|---|---|---|"]
    cidx = codes.set_index("taxonomy_code")
    for gid, g in cfg["groups"].items():
        for x in g.get("excluded_considered") or []:
            r = cidx.loc[x["code"]]
            lines.append(f"| {gid} | {x['code']} | {x['display_name']} | {int(r.n_npis_any):,} | "
                         f"{int(r.n_npis_primary):,} | {x['reason']} |")
    lines += ["", "Key ambiguity notes:", ""]
    for gid in ["primary_care", "electrophysiology", "neurology", "vascular", "pmr", "clinical_research"]:
        lines.append(f"* **{gid}** — {cfg['groups'][gid]['notes']}")
    lines += ["", "## Key variables (table `providers`)", "",
              "| column | meaning | used downstream for |", "|---|---|---|",
              "| npi | 10-digit National Provider Identifier (string) | key; `source_record_id` |",
              "| entity_type / entity_type_label | 1 individual, 2 organisation | individual vs facility counts |",
              "| name, provider_last_name, provider_first_name, organization_name, credential | legal names; credential "
              "is free text (individuals only) | display |",
              "| primary_taxonomy_code (+ NUCC display/grouping/classification/specialization) | taxonomy flagged "
              "primary | primary specialty |",
              "| all_taxonomy_codes | all distinct codes (pipe-separated) | audit |",
              "| specialty_groups / primary_specialty_groups / group_types | groups matched by any / primary code | "
              "clinic matching (`clinical_specialty_match`) |",
              "| address_line1, city, state, zip5, phone | primary practice location (self-reported) | display, "
              "geocoding |",
              "| lat, lon, county_fips, geocode_method | ZIP5-as-ZCTA internal point and containing 2024 county | "
              "maps, radius search, county density |",
              "| enumeration_date, last_update_date, deactivation_date, reactivation_date, certification_date | NPPES dates "
              "| record staleness |",
              "| n_secondary_practice_locations | rows in the practice-location reference file | flags multi-site providers |",
              "", "Other outputs: `provider_specialty_groups` (NPI × group, with primary flag and matched codes), "
              "`provider_density_county` (2024 county × group counts incl. zeros), `provider_taxonomy_counts` "
              "(national counts per group and per considered taxonomy code), `nucc_taxonomy`.", "",
              "## Missingness (measured on `providers`)", "",
              "| variable | denominator rows | missing (n, %) |", "|---|---|---|"]
    lines += _missingness(prov, ["npi", "name", "primary_taxonomy_code", "primary_taxonomy_display_name",
                                 "address_line1", "city", "state", "zip5", "phone", "enumeration_date",
                                 "last_update_date", "lat", "lon", "county_fips"])
    lines += [f"| credential (individuals only) | {int(ind.sum()):,} | "
              f"{_pct(int(prov.loc[ind, 'credential'].isna().sum()), int(ind.sum()))} |",
              f"| certification_date | {n:,} | {_pct(int(prov['certification_date'].isna().sum()), n)} |", "",
              f"* Geocoded (ZIP5 matches a 2024 ZCTA): {_pct(int(geocoded.sum()), n)}. Not geocoded by reason: "
              f"{prov.loc[~geocoded, 'not_geocoded_reason'].value_counts().to_dict()} (PO-box/unique ZIPs have no ZCTA).",
              "* Ungeocoded providers are concentrated in a few non-ZCTA ZIP5s (typically single-institution unique "
              "ZIPs), so county specialist counts are biased LOW in the counties that hold them. Top 10 (ZIP5: NPIs, "
              f"modal city/state): {_top_ungeocoded_zips(prov[~geocoded])}. Ungeocoded share by SPEC group: "
              + ", ".join(f"{r.specialty_group} {100.0 * r.n_not_geocoded / r.n_npis_any:.1f}%"
                          for r in grp.itertuples() if r.spec_listed) + ".",
              f"* County assignment among geocoded rows: "
              f"{prov.loc[geocoded, 'county_assignment_method'].value_counts(dropna=False).to_dict()}. "
              "`zcta_largest_land_share_rel2020` = the ZCTA's Census internal point lies outside every 2024 county "
              "polygon (coastal ZCTAs whose point falls in water, e.g. 60611 Chicago), so the county with the largest "
              "land share in the 2020 ZCTA-county relationship file is used (only if it is a 2024 county). "
              f"ZIPs affected: {fb_zips}.",
              f"* Geocoded rows whose assigned county lies in a different state than the address state: "
              f"{_pct(int((prov['geo_state_matches_address_state'] == False).sum()), int(geocoded.sum()))}.",  # noqa: E712
              "* **Connecticut vintage.** Counties are the 2024 vintage (CT = 9 planning regions 09110-09190); the "
              "relationship-file fallback only assigns 2024 counties. Measured: "
              f"{outs['connecticut']['ct_address_providers']:,} providers with a CT address; "
              f"{outs['connecticut']['ct_geocoded_to_2024_planning_region']:,} placed in a 2024 planning region, "
              f"{outs['connecticut']['ct_geocoded_to_legacy_county']:,} on a legacy CT county (09001-09015), "
              f"{outs['connecticut']['ct_geocoded_to_other_state_county']:,} in another state's county (ZIP/state "
              f"mismatch), {outs['connecticut']['ct_not_geocoded']:,} not geocoded. Any provider on a legacy CT "
              f"county: {outs['connecticut']['any_provider_on_legacy_ct_county']}. Tables on legacy CT counties "
              "(e.g. CMS) will not join `provider_density_county` for CT; report, do not force a mapping.",
              f"* Geocoded NPIs whose county is not in the 2024 county list: {len(outs['unmatched_counties'])} counties "
              f"({outs['unmatched_counties'][:10]}).",
              f"* Record staleness — last_update_date quantiles (10/25/50/75/90%) over all kept US NPIs: "
              f"{stats['kept_last_update_year_quantiles']}. Old dates mean the record has not been touched, "
              "not that it is wrong, but addresses may be out of date.",
              f"* Providers with ≥1 secondary practice location: "
              f"{_pct(int((prov['n_secondary_practice_locations'] > 0).sum()), n)}.", "",
              "## Linkage strategy", "",
              "* `npi` joins `providers` ↔ `provider_specialty_groups`; `county_fips` (2024 vintage) joins "
              "`provider_density_county` to `geographies` and to county-level geographic tables — an **ecological** join.",
              "* NPIs can be joined to other NPI-keyed public sources (e.g. CMS utilisation files) but NOT to any "
              "person-level participant data: providers are not participants, and nothing here says which patients a "
              "provider sees.",
              "* ClinicalTrials.gov / RePORTER sites carry no NPI; matching them to NPPES organisations would be fuzzy "
              "name/address matching and must be labelled as such.", "",
              "## Limitations and caveats", "",
              "* " + CAVEAT,
              "* Taxonomy codes are self-selected and need not be kept current (median last_update_date of kept US "
              f"NPIs: {stats['kept_last_update_year_quantiles'][2]}); a code shows self-description, not current "
              "practice, capacity, or acceptance of new patients.",
              "* Counts are NPIs, not FTEs: an individual clinician and the organisations they work for each hold "
              "NPIs, and organisations can register subparts (e.g. hospital departments) as separate NPIs "
              f"({_pct(int((prov['is_organization_subpart'] == 'Y').sum()), int((~ind).sum()))} of organisations in "
              "`providers` are subparts). A practice address can be an administrative/billing address.",
              "* The primary-care group includes 207R00000X Internal Medicine, which is also carried by hospital-based "
              f"and subspecialty internists, so it is an upper bound on outpatient primary care. {im_note} Nurse "
              "practitioners are included only for primary-care population foci and physician assistants are excluded "
              "(NUCC does not encode PA specialty).",
              "* Facility groups are self-reported NPPES taxonomies, not certification rosters (use HRSA for FQHC sites, "
              "CMS certification files for RHCs/CAHs).",
              "* ZIP-centroid geocoding can place a practice in the wrong county for ZIPs that cross county lines; "
              "PO-box/unique ZIPs are not geocoded.",
              f"* Territories are included; ZCTA coverage there is limited ({geo_rate}).", "",
              "## Processed outputs", "", "| table | rows |", "|---|---|",
              f"| providers | {len(outs['providers']):,} |",
              f"| provider_specialty_groups | {len(outs['membership']):,} |",
              f"| provider_density_county | {len(outs['density']):,} |",
              f"| provider_taxonomy_counts | {len(outs['counts']):,} |", "",
              "## Reproduce", "", "`uv run python -m measure_it.ingestion.nppes` (add `--refresh-listing` to pick up a "
              "newer monthly file). Stats: `data/raw/nppes/ingest_stats.json`; format comparison: "
              "`data/raw/nppes/format_version_comparison.json`."]
    p = raw_dir(SOURCE_ID) / "DATA_AUDIT.md"
    p.write_text("\n".join(lines) + "\n")
    return p


def write_nucc_audit(nucc: pd.DataFrame, rel: dict, cfg: dict, stats: dict) -> Path:
    e = load_manifest(NUCC_SOURCE_ID)["files"][rel["csv_file"]]
    sec = nucc["Section"].value_counts().to_dict()
    lines = ["# DATA AUDIT — NUCC Health Care Provider Taxonomy Code Set", "",
             "| Field | Value |", "|---|---|",
             f"| source_id | `{NUCC_SOURCE_ID}` |",
             f"| Source (dataset/API name, exact files/endpoints) | NUCC taxonomy CSV `{rel['csv_url']}` (linked as the "
             f"current version from `{rel['page_url']}`) |",
             "| Publishing organization | National Uniform Claim Committee (NUCC) |",
             f"| Retrieval date (UTC) | {e['retrieved_at']} |",
             f"| Source version / release | Version {rel['version']} (effective {rel.get('effective')}) |",
             f"| Source update date / cadence | Last-Modified {e.get('last_modified')}; NUCC releases twice a year "
             "(January and July) |",
             "| License / access conditions | Free public download from nucc.org; NUCC terms of use apply; "
             "no registration. |",
             "| Unit of observation | Taxonomy code |",
             f"| Sample size (actual, as ingested) | {len(nucc):,} codes ({sec}) |",
             "| Geography (resolution, vintage) | none |",
             "| Person-level? | no |", "| Geographic? | no |", "| Omics? | no |", "| Wearable? | no |",
             "| True participant linkage across modalities? | not applicable |", "",
             "## Files / endpoints retrieved", "",
             f"* `{rel['csv_file']}` — {e['bytes']:,} bytes, sha256 `{e['sha256']}` (MANIFEST.json). Previous versions "
             f"linked on the page: {', '.join(rel['previous_versions_linked'][:6])} …", "",
             "## Key variables", "",
             "Code (10 chars) → Grouping → Classification → Specialization; Display Name; Definition; Notes; Section "
             "(Individual / Non-Individual). Used to name NPPES taxonomy codes and to build "
             "`configs/specialty_groups.yaml`; `role` is derived (physician / nurse_practitioner / physician_assistant "
             "/ other_individual / facility).", "",
             "## Missingness", "", "| column | missing (n, %) |", "|---|---|"]
    for c in nucc.columns:
        miss = int((nucc[c].astype(str).str.strip() == "").sum())
        lines.append(f"| {c} | {_pct(miss, len(nucc))} |")
    ndef = int(nucc["Definition"].str.contains("Definition to come", case=False).sum())
    lines += [f"| Definition = 'Definition to come...' | {_pct(ndef, len(nucc))} |", "",
              "## Linkage strategy", "",
              f"Joins to NPPES on the taxonomy code. Coverage: of {stats['distinct_taxonomy_codes']:,} distinct codes used "
              f"by active US NPPES records, codes absent from this NUCC version: {stats['codes_not_in_nucc']}.", "",
              "## Limitations and caveats", "",
              "* A taxonomy is a self-selected provider classification, not a credential, board certification, or "
              "statement of which conditions a provider treats.",
              "* Many facility codes have no definition ('Definition to come...').",
              "* Deactivated codes are not in the current CSV; NPPES records may still carry them (see coverage above).", "",
              "## Processed outputs", "", f"| nucc_taxonomy | {len(nucc):,} |", "",
              "## Reproduce", "", "`uv run python -m measure_it.ingestion.nppes` (writes both sources)."]
    p = raw_dir(NUCC_SOURCE_ID) / "DATA_AUDIT.md"
    p.write_text("\n".join(lines) + "\n")
    return p


def write_registry(stats: dict, outs: dict, release: dict, meta: dict, nucc_rel: dict, nucc_rows: int) -> None:
    grp = outs["counts"]
    grp = grp[grp["count_level"] == "specialty_group"].set_index("specialty_group")["n_npis_any"]
    write_registry_entry({
        "source_id": SOURCE_ID, "name": "NPPES NPI Registry downloadable file (V.2, full replacement monthly)",
        "publisher": "CMS (NPPES)", "landing_url": LISTING_URL,
        "access_urls": [release["monthly_v2_url"], release["deactivation_url"]],
        "license": "FOIA-disclosable U.S. government data (public)",
        "access_conditions": "open download, no registration or data-use agreement",
        "retrieved_at": meta["retrieved_at"], "source_version": meta["source_version"],
        "update_date": meta.get("last_modified") or "monthly",
        "data_layer": "facility", "unit_of_observation": "NPI (individual provider or organisation)",
        "sample_size": {"npi_rows_in_file": int(stats["rows_total"]),
                        "deactivated_excluded": int(stats["deactivated_excluded"]),
                        "active_us_practice_location": int(stats["kept_us_rows"]),
                        "providers_in_any_group": int(len(outs["providers"])),
                        "per_group_npis_any": {k: int(v) for k, v in grp.items()}},
        "geographic_resolution": "practice ZIP5 -> ZCTA internal point (zcta_centroid); 2024 county aggregates",
        "person_level": False, "geographic": True, "omics": False, "wearable": False,
        "participant_linkage": "not applicable (providers, not participants)",
        "true_participant_linkage_across_modalities": False, "status": "ingested",
        "processed_outputs": ["providers", "provider_specialty_groups", "provider_density_county",
                              "provider_taxonomy_counts"],
        "audit": f"data/raw/{SOURCE_ID}/DATA_AUDIT.md", "ingestion_module": "measure_it.ingestion.nppes",
        "limitations": [
            "NPI taxonomy is self-reported and does not mean a provider treats Long COVID, ME/CFS or POTS",
            "NPPES addresses are self-reported and may be billing-oriented",
            "geocodes are ZIP centroids (ZCTA internal points), not address geocodes",
            "counts are NPIs, not FTEs; organisations register subparts separately",
            "primary practice location only; secondary locations counted, not mapped",
        ],
        "format_version": "V.2 (V.1 retired by CMS 2026-03-03; header columns identical, longer name fields)",
    })
    write_registry_entry({
        "source_id": NUCC_SOURCE_ID, "name": "NUCC Health Care Provider Taxonomy code set (CSV)",
        "publisher": "National Uniform Claim Committee", "landing_url": nucc_rel["page_url"],
        "access_urls": [nucc_rel["csv_url"]], "license": "free public download; NUCC terms of use",
        "access_conditions": "open download, no registration",
        "retrieved_at": load_manifest(NUCC_SOURCE_ID)["files"][nucc_rel["csv_file"]]["retrieved_at"],
        "source_version": nucc_rel["version"], "update_date": nucc_rel.get("effective") or "semi-annual",
        "data_layer": "facility", "unit_of_observation": "taxonomy code", "sample_size": {"codes": int(nucc_rows)},
        "geographic_resolution": "none", "person_level": False, "geographic": False, "omics": False,
        "wearable": False, "participant_linkage": "not applicable",
        "true_participant_linkage_across_modalities": False, "status": "ingested",
        "processed_outputs": ["nucc_taxonomy"], "audit": f"data/raw/{NUCC_SOURCE_ID}/DATA_AUDIT.md",
        "ingestion_module": "measure_it.ingestion.nppes",
        "limitations": ["a taxonomy is a self-selected classification, not a credential or a list of conditions treated",
                        "many facility codes have no NUCC definition"],
    })


# ---------------------------------------------------------------- entry point
def run(refresh_listing: bool = False, skip_format_compare: bool = False) -> dict:
    release = discover_nppes_release(refresh=refresh_listing)
    (raw_dir(SOURCE_ID) / "release.json").write_text(json.dumps(release, indent=2))
    zpath = download_file(release["monthly_v2_url"], SOURCE_ID, release["monthly_v2_file"],
                          timeout=1800, min_bytes=100_000_000)
    deact_path = download_file(release["deactivation_url"], SOURCE_ID, release["deactivation_file"]) \
        if release["deactivation_url"] else None
    zentry = load_manifest(SOURCE_ID)["files"][release["monthly_v2_file"]]

    nucc, nucc_rel = _nucc_frame()
    (raw_dir(NUCC_SOURCE_ID) / "release.json").write_text(json.dumps(nucc_rel, indent=2))
    cfg = sp.load_specialty_groups()
    problems = sp.validate_against_nucc(cfg, nucc)
    if problems:
        raise ValueError("specialty_groups.yaml does not match NUCC file: " + "; ".join(problems))
    nucc_retrieved = load_manifest(NUCC_SOURCE_ID)["files"][nucc_rel["csv_file"]]["retrieved_at"]
    write_table(build_nucc_table(nucc, cfg, nucc_rel, nucc_retrieved), "nucc_taxonomy", producer=PRODUCER,
                description=f"NUCC provider taxonomy v{nucc_rel['version']} with specialty-group membership")

    members = extract_members(zpath, raw_dir(SOURCE_ID) / "extracted")
    meta = {
        "retrieved_at": zentry["retrieved_at"], "last_modified": zentry.get("last_modified"),
        "csv_name": members["npidata"].name, "csv_bytes": members["npidata"].stat().st_size,
        "pl_name": members["pl"].name, "nucc_version": nucc_rel["version"],
        "source_version": f"{release['monthly_v2_file']} ({members['npidata'].name}); NUCC {nucc_rel['version']}",
    }
    stats, frames, con, csv_sql = process(members, cfg, nucc)
    deact = deactivation_crosscheck(deact_path, con) if deact_path else {}
    fmt = {} if skip_format_compare else compare_format_versions(members, con, csv_sql)
    con.close()

    outs = build_outputs(frames, cfg, nucc, meta)
    write_table(outs["providers"], "providers", producer=PRODUCER,
                description="NPPES NPIs (active, US practice location) in >=1 specialty/facility group")
    write_table(outs["membership"], "provider_specialty_groups", producer=PRODUCER,
                description="NPI x specialty/facility group membership (any/primary taxonomy)")
    write_table(outs["density"], "provider_density_county", producer=PRODUCER,
                description="2024 county x group NPI counts (no denominators; zeros kept)")
    write_table(outs["counts"], "provider_taxonomy_counts", producer=PRODUCER,
                description="National NPI counts per group and per considered taxonomy code")

    stats["deactivation_crosscheck"] = deact
    stats["unmatched_counties"] = outs["unmatched_counties"]
    stats["connecticut"] = outs["connecticut"]
    stats["output_rows"] = {k: int(len(outs[k])) for k in ("providers", "membership", "density", "counts")}
    (raw_dir(SOURCE_ID) / "ingest_stats.json").write_text(json.dumps(stats, indent=2, default=str))
    write_nppes_audit(stats, outs, cfg, release, fmt, deact, zentry, meta)
    write_nucc_audit(nucc, nucc_rel, cfg, stats)
    write_registry(stats, outs, release, meta, nucc_rel, len(nucc))
    return stats


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description="Ingest NPPES (V.2) + NUCC taxonomy")
    ap.add_argument("--refresh-listing", action="store_true", help="re-read the CMS/NUCC listing pages")
    ap.add_argument("--skip-format-compare", action="store_true", help="skip the V.1 vs V.2 layout comparison")
    a = ap.parse_args(argv)
    stats = run(refresh_listing=a.refresh_listing, skip_format_compare=a.skip_format_compare)
    print(json.dumps({k: stats[k] for k in ("rows_total", "deactivated_excluded", "kept_us_rows", "output_rows")},
                     indent=2, default=str))


if __name__ == "__main__":
    main()
