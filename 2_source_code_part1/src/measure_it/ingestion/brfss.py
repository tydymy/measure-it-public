"""CDC BRFSS annual LLCP public files: Long COVID items, design variables and covariates (2022-2024).

Source: Behavioral Risk Factor Surveillance System, combined landline + cell-phone ("LLCP") public data,
SAS transport files from https://www.cdc.gov/brfss/annual_data/annual_<year>.html, plus each year's codebook
(USCODE<yy>_LLCP HTML) and one-column variable layout, and the SMART BRFSS MMSA file for the same year.

What the files hold (verified from the codebooks/layouts on every run, never assumed):
  2022  Long-term COVID Effects core section: COVIDPOS (told tested positive; 3 = home test only),
        COVIDSMP ("Did you have any symptoms lasting 3 months or longer ..." = EVER long COVID), COVIDPRM.
  2023  Long-term COVID Effects core section: COVIDPO1, COVIDSM1 ("Do you currently have symptoms lasting 3
        months or longer ..." = CURRENT long COVID), COVIDACT (activity limitation: a lot / a little / not at all).
        Kentucky and Pennsylvania are absent from the 2023 public file (CDC: insufficient data).
  2024  no COVID item in the layout (checked).
Design: _STSTR, _PSU, _LLCPWT. Geography on the public LLCP file: _STATE and two county-derived urbanicity
flags (_METSTAT, _URBSTAT, from the NCHS urban-rural code of the respondent's county); no county identifier.
The SMART MMSA file adds _MMSA / _MMSAWT for metropolitan areas with >= 500 respondents.

Validation: 2023 direct state estimates of current Long COVID (Taylor-linearised, strata x PSU) are compared
with the MMWR supplementary table (Ford et al., MMWR 2024;73(50), PMC11658399; stacks.cdc.gov/view/cdc/174567).

Outputs (data/processed): brfss_long_covid_direct_estimates (national + state direct survey estimates; not a
geo_condition_burden partition: the coordinator decides whether state BRFSS rows join the burden table).
Person-level microdata stay in data/interim/cdc_brfss (never written as a processed person table, never joined to
geography beyond the survey's own state/urbanicity fields).

Run: uv run python -m measure_it.ingestion.brfss
"""
from __future__ import annotations

import html
import json
import re
import shutil
import subprocess
import zipfile
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd

from ..config import RAW, RESULTS, interim_dir, raw_dir
from ..download import download_file, load_manifest
from ..provenance import add_provenance
from ..registry import write_registry_entry
from ..store import read_table, write_table

SOURCE_ID = "cdc_brfss"
SOURCE_NAME = "CDC Behavioral Risk Factor Surveillance System (BRFSS) LLCP public data"
PRODUCER = "measure_it.ingestion.brfss"
YEARS = (2022, 2023, 2024)
BASE = "https://www.cdc.gov/brfss/annual_data/{y}"
LANDING = "https://www.cdc.gov/brfss/annual_data/annual_{y}.html"
SMART_LANDING = "https://www.cdc.gov/brfss/smart/smart_{y}.html"
MMWR_URL = "https://pmc.ncbi.nlm.nih.gov/articles/PMC11658399/"
MMWR_SUPP_URL = "https://stacks.cdc.gov/view/cdc/174567/cdc_174567_DS1.pdf"
MMWR_SUPP_FILE = "validation/mmwr_7350a2_supplementary_table_cdc_174567_DS1.pdf"
NCHS_URC_URL = "https://www.cdc.gov/nchs/data/data_acces_files/NCHSURCodes2013.xlsx"
OMB_2020_URL = ("https://www2.census.gov/programs-surveys/metro-micro/geographies/reference-files/2020/"
                "delineation-files/list1_2020.xls")
OMB_2023_URL = ("https://www2.census.gov/programs-surveys/metro-micro/geographies/reference-files/2023/"
                "delineation-files/list1_2023.xlsx")
ACS_SF = "https://www2.census.gov/programs-surveys/acs/summary_file/2024/table-based-SF/data/5YRData"
ACS_POSTSTRAT_TABLES = ("B01001B", "B01001H", "B01001I", "B03002", "B15001")

# Long-COVID items expected per year (checked against layout + codebook; an absent item raises, never guessed).
LC_VARS = {2022: ["COVIDPOS", "COVIDSMP", "COVIDPRM"], 2023: ["COVIDPO1", "COVIDSM1", "COVIDACT"], 2024: []}
DESIGN_VARS = ["_STATE", "_STSTR", "_PSU", "_LLCPWT"]
COVARIATES = ["_AGEG5YR", "_SEX", "_IMPRACE", "_EDUCAG", "_INCOMG1", "_URBSTAT", "_METSTAT", "DISPCODE", "IYEAR"]
MMSA_VARS = ["_MMSA", "_MMSAWT", "MMSANAME", "_STSTR", "_AGEG5YR", "_SEX", "DISPCODE"]
SUBSTATE_PATTERN = re.compile(r"CNTY|COUNTY|CTYCODE|_MMSA|ZIP|TRACT", re.IGNORECASE)

TERRITORY_FIPS = {"Guam": "66", "Puerto Rico": "72", "Virgin Islands": "78"}
STATES_50_DC = None  # resolved lazily from `geographies`

# Model / post-stratification categories (shared with geography.sae).
AGE7 = ["18-24", "25-34", "35-44", "45-54", "55-64", "65-74", "75+"]
AGEG5YR_TO_AGE7 = {1: 0, 2: 1, 3: 1, 4: 2, 5: 2, 6: 3, 7: 3, 8: 4, 9: 4, 10: 5, 11: 5, 12: 6, 13: 6}
AGE3 = ["18-44", "45-64", "65+"]
AGEG5YR_TO_AGE3 = {**{k: 0 for k in range(1, 6)}, **{k: 1 for k in range(6, 10)}, **{k: 2 for k in range(10, 14)}}
SEX2 = ["male", "female"]
RACE4 = ["white_nh", "black_nh", "hispanic", "other_nh"]
IMPRACE_TO_RACE4 = {1: 0, 2: 1, 5: 2, 3: 3, 4: 3, 6: 3}
EDU4 = ["lt_hs", "hs", "some_college", "college_grad"]
URBAN3 = ["metro", "micropolitan", "noncore"]

MEASURES = {
    # measure_id: (year, label, denominator description)
    "lc_current_pct_all_adults": (2023, "Current long COVID (COVIDPO1=1 and COVIDSM1=1), % of adults answering",
                                  "adults 18+ with a valid answer (complete case: partial interviews that stopped "
                                  "before the COVID section, don't know and refused excluded)"),
    "lc_current_pct_all_adults_mmwr_def": (2023, "Current long COVID, MMWR definition (non-response counted as no)",
                                           "all adults 18+ in the file (MMWR 2024;73(50) denominator)"),
    "lc_significant_limitation_pct_all_adults": (2023, "Current long COVID with activity reduced 'a lot' (COVIDACT=1)",
                                                 "adults 18+ with a valid answer (complete case)"),
    "lc_significant_limitation_pct_current_lc_mmwr_def": (2023, "Significant activity limitation among adults "
                                                         "with current long COVID (MMWR definition)",
                                                         "adults with current long COVID (COVIDSM1=1)"),
    "lc_ever_pct_all_adults": (2022, "Ever had symptoms >= 3 months after COVID (COVIDPOS in 1,3 and COVIDSMP=1)",
                               "adults 18+ with a valid answer (complete case)"),
    "lc_ever_pct_all_adults_all_records": (2022, "Ever long COVID, non-response counted as no",
                                           "all adults 18+ in the file"),
}


# --------------------------------------------------------------------------- urls
def year_urls(y: int) -> dict[str, tuple[str, str]]:
    yy = str(y)[2:]
    b = BASE.format(y=y)
    return {
        "xpt": (f"{b}/files/LLCP{y}XPT.zip", f"{y}/LLCP{y}XPT.zip"),
        "codebook": (f"{b}/zip/codebook{yy}_llcp-v2-508.zip", f"{y}/codebook{yy}_llcp-v2-508.zip"),
        "layout": (f"{b}/llcp_varlayout_{yy}_onecolumn.html", f"{y}/llcp_varlayout_{yy}_onecolumn.html"),
        "overview": (f"{b}/pdf/Overview_{y}-508.pdf", f"{y}/Overview_{y}-508.pdf"),
        "mmsa_layout": (f"{b}/pdf/MMSA_VarLayout_{yy}-508.pdf", f"{y}/mmsa/MMSA_VarLayout_{yy}-508.pdf"),
        "mmsa_method": (f"{b}/pdf/{y}_SMART_BRFSS_MMSA_Methodology-508.pdf",
                        f"{y}/mmsa/{y}_SMART_BRFSS_MMSA_Methodology-508.pdf"),
        "mmsa_xpt": (f"{b}/files/MMSA{y}_XPT.zip", f"{y}/mmsa/MMSA{y}_XPT.zip"),
    }


def fetch(years=YEARS, mmsa_years=(2023,)) -> dict[str, Path]:
    out = {}
    for y in years:
        for key, (url, fn) in year_urls(y).items():
            if key.startswith("mmsa") and y not in mmsa_years:
                continue
            out[f"{y}:{key}"] = download_file(url, SOURCE_ID, fn, allow_html=fn.endswith(".html"), timeout=1800,
                                              max_retries=3, min_bytes=1000)
    out["mmwr_supp"] = download_file(MMWR_SUPP_URL, SOURCE_ID, MMWR_SUPP_FILE)
    out["nchs_urc"] = download_file(NCHS_URC_URL, SOURCE_ID, "reference/NCHSURCodes2013.xlsx")
    out["omb_2020"] = download_file(OMB_2020_URL, SOURCE_ID, "reference/omb_delineation_list1_2020.xls")
    out["omb_2023"] = download_file(OMB_2023_URL, SOURCE_ID, "reference/omb_delineation_list1_2023.xlsx")
    for t in ACS_POSTSTRAT_TABLES:
        out[f"acs:{t}"] = download_file(f"{ACS_SF}/acsdt5y2024-{t.lower()}.dat", SOURCE_ID,
                                        f"acs_poststrat/acsdt5y2024-{t.lower()}.dat", min_bytes=10_000,
                                        timeout=1800, max_retries=3)
    return out


# --------------------------------------------------------------------------- codebook / layout parsing (pure)
def _strip_html(text: str) -> str:
    t = re.sub(r"<[^>]+>", " ", text)
    t = html.unescape(t).replace("\xa0", " ")
    return re.sub(r"\s+", " ", t)


def layout_variables(layout_html: str) -> list[str]:
    """Variable names from the one-column layout page (table cells that look like SAS names)."""
    cells = [_strip_html(c).strip() for c in re.findall(r"<td[^>]*>(.*?)</td>", layout_html, re.S)]
    return [c for c in cells if re.fullmatch(r"_?[A-Z][A-Z0-9_]+", c)]


def codebook_entries(codebook_html: str) -> dict[str, dict]:
    """SAS variable name -> {label, section, question, text} from the USCODE codebook HTML."""
    t = _strip_html(codebook_html)
    out = {}
    for part in re.split(r"(?= Label: )", t):
        m = re.search(r"SAS Variable Name: (\S+)", part)
        if not m:
            continue
        label = re.search(r"Label: (.*?) Section Name:", part)
        section = re.search(r"Section Name: (.*?) (?:Section|Module) Number:", part)
        question = re.search(r"Question: (.*?) Value Value Label", part)
        out[m.group(1)] = {"label": label.group(1).strip() if label else "",
                           "section": section.group(1).strip() if section else "",
                           "question": question.group(1).strip() if question else "",
                           "text": part.strip()}
    return out


def read_codebook(zip_path: Path) -> tuple[str, dict[str, dict]]:
    with zipfile.ZipFile(zip_path) as z:
        names = [n for n in z.namelist() if n.lower().endswith((".html", ".htm"))]
        raw = z.read(names[0]).decode("latin-1")
    return names[0], codebook_entries(raw)


# --------------------------------------------------------------------------- microdata
def read_subset(zip_path: Path, columns: list[str], cache: Path) -> pd.DataFrame:
    """Read selected columns of the XPT inside a zip (extract transiently to interim, cache as Parquet)."""
    if cache.exists():
        return pd.read_parquet(cache)
    import pyreadstat
    with zipfile.ZipFile(zip_path) as z:
        member = [n for n in z.namelist() if n.strip().lower().endswith(".xpt")][0]
        tmp = cache.with_suffix(".xpt.part")
        with z.open(member) as src, open(tmp, "wb") as dst:
            shutil.copyfileobj(src, dst, 1 << 24)
    try:
        _, meta = pyreadstat.read_xport(str(tmp), metadataonly=True, encoding="latin1")
        cols = [c for c in columns if c in meta.column_names]
        df, _ = pyreadstat.read_xport(str(tmp), usecols=cols, encoding="latin1")
    finally:
        tmp.unlink(missing_ok=True)
    df.to_parquet(cache, index=False)
    return df


def _code(s: pd.Series, mapping: dict) -> pd.Series:
    return s.map(mapping).astype("Float64")


def derive(df: pd.DataFrame, year: int) -> pd.DataFrame:
    """Analysis variables (outcomes as 1/0/NaN; categories as integer codes, NaN = not post-stratifiable)."""
    out = pd.DataFrame(index=df.index)
    out["state_fips"] = df["_STATE"].astype(int).astype(str).str.zfill(2)
    out["ststr"] = df["_STSTR"].astype("int64")
    out["psu"] = df["_PSU"].astype("int64") if "_PSU" in df else np.arange(len(df))
    out["weight"] = df["_LLCPWT"].astype(float)
    out["year"] = year
    out["partial_interview"] = df["DISPCODE"].eq(1200)
    out["age7"] = _code(df["_AGEG5YR"], AGEG5YR_TO_AGE7)
    out["age3"] = _code(df["_AGEG5YR"], AGEG5YR_TO_AGE3)
    out["sex"] = _code(df["_SEX"], {1: 0, 2: 1})
    out["race4"] = _code(df["_IMPRACE"], IMPRACE_TO_RACE4)
    out["edu4"] = _code(df["_EDUCAG"], {1: 0, 2: 1, 3: 2, 4: 3})
    met, urb = df["_METSTAT"], df["_URBSTAT"]
    out["urban3"] = pd.Series(np.select([met.eq(1), met.eq(2) & urb.eq(1), urb.eq(2)], [0, 1, 2], default=-1),
                              index=df.index).replace(-1, np.nan).astype("Float64")
    nan = np.nan
    if year == 2023:
        po, sm, act = df["COVIDPO1"], df["COVIDSM1"], df["COVIDACT"]
        case = po.eq(1) & sm.eq(1)
        valid_no = po.eq(2) | (po.eq(1) & sm.eq(2))
        out["lc_current_pct_all_adults"] = np.where(case, 1.0, np.where(valid_no, 0.0, nan))
        out["lc_current_pct_all_adults_mmwr_def"] = case.astype(float)
        sig_known = case & act.isin([1, 2, 3])
        out["lc_significant_limitation_pct_all_adults"] = np.where(
            case & act.eq(1), 1.0, np.where(valid_no | sig_known, 0.0, nan))
        out["lc_significant_limitation_pct_current_lc_mmwr_def"] = np.where(case, act.eq(1).astype(float), nan)
    elif year == 2022:
        pos, smp = df["COVIDPOS"], df["COVIDSMP"]
        case = pos.isin([1, 3]) & smp.eq(1)
        valid_no = pos.eq(2) | (pos.isin([1, 3]) & smp.eq(2))
        out["lc_ever_pct_all_adults"] = np.where(case, 1.0, np.where(valid_no, 0.0, nan))
        out["lc_ever_pct_all_adults_all_records"] = case.astype(float)
    return out


# --------------------------------------------------------------------------- design-based estimation (pure)
def weighted_prop(y, w, strata, psu, domain=None) -> dict:
    """Weighted proportion with Taylor-linearised SE (with-replacement PSUs within strata).

    Records with y = NaN or outside `domain` contribute zero to the linearised variable (domain estimation keeps
    the full design). Singleton strata are centred on the grand mean of PSU totals (R survey 'adjust').
    """
    y = np.asarray(y, float)
    w = np.asarray(w, float)
    d = ~np.isnan(y) if domain is None else (~np.isnan(y) & np.asarray(domain, bool))
    wd = np.where(d, w, 0.0)
    W = wd.sum()
    if W <= 0:
        return {"p": np.nan, "se": np.nan, "n": 0, "cases": 0}
    yy = np.where(d, y, 0.0)
    p = (wd * yy).sum() / W
    z = wd * (yy - p) / W
    t = pd.DataFrame({"h": np.asarray(strata), "j": np.asarray(psu), "z": z}).groupby(["h", "j"], sort=False)["z"].sum()
    t = t.reset_index()
    nh = t.groupby("h")["z"].transform("size")
    mh = t.groupby("h")["z"].transform("mean")
    grand = t["z"].mean()
    dev = np.where(nh > 1, t["z"] - mh, t["z"] - grand)
    fac = np.where(nh > 1, nh / (nh - 1).clip(lower=1), 1.0)
    var = float((fac * dev ** 2).sum())
    return {"p": float(p), "se": float(np.sqrt(var)), "n": int(d.sum()), "cases": int((yy[d] == 1).sum())}


def logit_ci(p: float, se: float, z: float = 1.959964) -> tuple[float, float]:
    """95% CI on the logit scale (as SUDAAN/survey 'logit' intervals)."""
    if not (0 < p < 1) or not np.isfinite(se):
        return (np.nan, np.nan)
    lse = se / (p * (1 - p))
    lo, hi = np.log(p / (1 - p)) - z * lse, np.log(p / (1 - p)) + z * lse
    return (float(1 / (1 + np.exp(-lo))), float(1 / (1 + np.exp(-hi))))


def standardized_prop(df: pd.DataFrame, ycol: str, standard: dict[tuple[int, int], float], domain=None) -> dict:
    """Direct age3 x sex standardised proportion; variance = sum s_g^2 var_g (cell covariances ignored)."""
    tot = sum(standard.values())
    p = v = 0.0
    for (a, s), share in standard.items():
        cell = (df["age3"] == a) & (df["sex"] == s)
        if domain is not None:
            cell &= domain
        r = weighted_prop(df[ycol], df["weight"], df["ststr"], df["psu"], domain=cell.fillna(False).to_numpy())
        if not np.isfinite(r["p"]):
            return {"p": np.nan, "se": np.nan}
        p += share / tot * r["p"]
        v += (share / tot) ** 2 * r["se"] ** 2
    return {"p": p, "se": float(np.sqrt(v))}


def acs_national_standard() -> dict[tuple[int, int], float]:
    """Adults by age3 x sex, United States, ACS 2020-2024 B01001 (proxy for MMWR's 2020 Census standard)."""
    path = RAW / "census_acs" / "acsdt5y2024-b01001.dat"
    con = duckdb.connect()
    row = con.execute("SELECT * FROM read_csv(?, delim='|', header=true, all_varchar=true) WHERE GEO_ID='0100000US'",
                      [str(path)]).df().iloc[0]
    con.close()
    # B01001 cells: male 007..025 (18-19 .. 85+), female 031..049; bounds per cell
    lows = [18, 20, 21, 22, 25, 30, 35, 40, 45, 50, 55, 60, 62, 65, 67, 70, 75, 80, 85]
    out = {}
    for s, start in ((0, 7), (1, 31)):
        for i, lo in enumerate(lows):
            a = 0 if lo < 45 else (1 if lo < 65 else 2)
            out[(a, s)] = out.get((a, s), 0.0) + float(row[f"B01001_E{start + i:03d}"])
    return out


# --------------------------------------------------------------------------- MMWR validation table (pure parser)
MMWR_ROW = re.compile(
    r"^\s*(?P<name>[A-Z][A-Za-z .]+?)\s+(?P<n1>[\d,]+)\s+(?P<p1>\d+\.\d)\s*\((?P<l1>\d+\.\d),\s*(?P<h1>\d+\.\d)\)"
    r"\s+(?P<n2>[\d,]+)\s+(?P<p2>\d+\.\d)\s*\((?P<l2>\d+\.\d),\s*(?P<h2>\d+\.\d)\)")


def parse_mmwr_table(text: str) -> pd.DataFrame:
    rows = []
    for line in text.splitlines():
        m = MMWR_ROW.match(line)
        if m:
            g = m.groupdict()
            rows.append({"jurisdiction": g["name"].strip(), "n_cases_current": int(g["n1"].replace(",", "")),
                         "mmwr_current_pct": float(g["p1"]), "mmwr_current_lo": float(g["l1"]),
                         "mmwr_current_hi": float(g["h1"]), "n_cases_sig": int(g["n2"].replace(",", "")),
                         "mmwr_sig_among_lc_pct": float(g["p2"]), "mmwr_sig_lo": float(g["l2"]),
                         "mmwr_sig_hi": float(g["h2"])})
        elif re.match(r"^\s*(Kentucky|Pennsylvania)\s+—", line):
            rows.append({"jurisdiction": line.split("—")[0].strip()})
    return pd.DataFrame(rows)


def pdf_text(path: Path) -> str | None:
    exe = shutil.which("pdftotext")
    if not exe:
        return None
    r = subprocess.run([exe, "-layout", "-q", str(path), "-"], capture_output=True, timeout=120, check=False)
    return r.stdout.decode("utf-8", errors="replace")


# --------------------------------------------------------------------------- estimates
def state_names() -> dict[str, str]:
    g = read_table("geographies")
    st = g[g["geo_level"] == "state"]
    return dict(zip(st["geo_id"], st["name"]))


def direct_estimates(d: pd.DataFrame, measure_ids: list[str], standard: dict) -> pd.DataFrame:
    rows = []
    names = state_names()
    for mid in measure_ids:
        groups = [("national", "US50DC", "United States (50 states + DC in file)", d["state_fips"].isin(
            [f for f in names if int(f) <= 56]))]
        groups += [("state", s, names.get(s, s), d["state_fips"].eq(s)) for s in sorted(d["state_fips"].unique())]
        for level, gid, gname, dom in groups:
            base_dom = dom.to_numpy()
            if mid == "lc_significant_limitation_pct_current_lc_mmwr_def":
                base_dom = base_dom & d[mid].notna().to_numpy()
            r = weighted_prop(d[mid], d["weight"], d["ststr"], d["psu"], domain=base_dom)
            lo, hi = logit_ci(r["p"], r["se"])
            rows.append({"measure_id": mid, "geo_level": level, "geo_id": gid, "geo_name": gname,
                         "estimate_type": "crude", "value": 100 * r["p"], "se": 100 * r["se"],
                         "ci_low": 100 * lo, "ci_high": 100 * hi, "n_unweighted": r["n"], "n_cases": r["cases"]})
            if mid.endswith("mmwr_def"):
                s = standardized_prop(d, mid, standard, domain=pd.Series(base_dom, index=d.index))
                lo, hi = logit_ci(s["p"], s["se"])
                rows.append({"measure_id": mid, "geo_level": level, "geo_id": gid, "geo_name": gname,
                             "estimate_type": "age_sex_standardized_acs2024", "value": 100 * s["p"],
                             "se": 100 * s["se"], "ci_low": 100 * lo, "ci_high": 100 * hi,
                             "n_unweighted": r["n"], "n_cases": r["cases"]})
    return pd.DataFrame(rows)


def mmwr_comparison(direct: pd.DataFrame, mmwr: pd.DataFrame) -> pd.DataFrame:
    names = state_names()
    name_to_fips = {v: k for k, v in names.items()} | TERRITORY_FIPS | {"National": "US50DC"}
    m = mmwr.copy()
    m["geo_id"] = m["jurisdiction"].map(name_to_fips)
    out = m
    for mid, col in (("lc_current_pct_all_adults_mmwr_def", "current"),
                     ("lc_significant_limitation_pct_current_lc_mmwr_def", "sig_among_lc")):
        for et, tag in (("crude", "crude"), ("age_sex_standardized_acs2024", "std")):
            sub = direct[(direct.measure_id == mid) & (direct.estimate_type == et)].set_index("geo_id")
            out[f"ours_{col}_{tag}"] = out["geo_id"].map(sub["value"])
            if tag == "std":
                out[f"ours_{col}_{tag}_lo"] = out["geo_id"].map(sub["ci_low"])
                out[f"ours_{col}_{tag}_hi"] = out["geo_id"].map(sub["ci_high"])
            out[f"ours_{col}_n_cases"] = out["geo_id"].map(sub["n_cases"])
    cc = direct[(direct.measure_id == "lc_current_pct_all_adults") & (direct.estimate_type == "crude")].set_index("geo_id")
    out["ours_current_complete_case_crude"] = out["geo_id"].map(cc["value"])
    return out


def agreement_summary(cmp_: pd.DataFrame) -> dict:
    s = cmp_.dropna(subset=["mmwr_current_pct", "ours_current_std"])
    s = s[s["geo_id"] != "US50DC"]
    diff = s["ours_current_std"] - s["mmwr_current_pct"]
    sig = cmp_.dropna(subset=["mmwr_sig_among_lc_pct", "ours_sig_among_lc_std"])
    sig = sig[sig["geo_id"] != "US50DC"]
    dsig = sig["ours_sig_among_lc_std"] - sig["mmwr_sig_among_lc_pct"]
    nat = cmp_[cmp_["geo_id"] == "US50DC"]
    cases_match = int((s["n_cases_current"] == s["ours_current_n_cases"]).sum())
    return {
        "n_jurisdictions": int(len(s)),
        "case_counts_identical": cases_match,
        "current_mae_pp": float(diff.abs().mean()), "current_max_abs_pp": float(diff.abs().max()),
        "current_within_0_1pp": int((diff.abs() <= 0.1 + 1e-9).sum()),
        "current_within_0_3pp": int((diff.abs() <= 0.3 + 1e-9).sum()),
        "current_pearson": float(np.corrcoef(s["ours_current_std"], s["mmwr_current_pct"])[0, 1]),
        "current_mae_crude_pp": float((s["ours_current_crude"] - s["mmwr_current_pct"]).abs().mean()),
        "sig_mae_pp": float(dsig.abs().mean()), "sig_max_abs_pp": float(dsig.abs().max()),
        "sig_pearson": float(np.corrcoef(sig["ours_sig_among_lc_std"], sig["mmwr_sig_among_lc_pct"])[0, 1]),
        "national_ours_std": float(nat["ours_current_std"].iloc[0]) if len(nat) else np.nan,
        "national_ours_crude": float(nat["ours_current_crude"].iloc[0]) if len(nat) else np.nan,
        "national_mmwr": float(nat["mmwr_current_pct"].iloc[0]) if len(nat) else np.nan,
        "national_cases_mmwr": int(nat["n_cases_current"].iloc[0]) if len(nat) else -1,
        "national_cases_ours": int(nat["ours_current_n_cases"].iloc[0]) if len(nat) else -1,
        "worst": s.assign(diff=diff).reindex(diff.abs().sort_values(ascending=False).index)[
            ["jurisdiction", "mmwr_current_pct", "ours_current_std", "diff"]].head(5).to_dict("records"),
    }


# --------------------------------------------------------------------------- MMSA
def mmsa_subset(year: int = 2023) -> pd.DataFrame:
    url, fn = year_urls(year)["mmsa_xpt"]
    zp = download_file(url, SOURCE_ID, fn, timeout=1800, max_retries=3)
    cols = MMSA_VARS + LC_VARS[year]
    return read_subset(zp, cols, interim_dir(SOURCE_ID) / f"mmsa{year}_subset.parquet")


def mmsa_direct(year: int = 2023) -> pd.DataFrame:
    m = mmsa_subset(year)
    po, sm = m["COVIDPO1"], m["COVIDSM1"]
    case = po.eq(1) & sm.eq(1)
    valid_no = po.eq(2) | (po.eq(1) & sm.eq(2))
    y = np.where(case, 1.0, np.where(valid_no, 0.0, np.nan))
    rows = []
    for code, g in m.assign(y=y).groupby("_MMSA"):
        r = weighted_prop(g["y"], g["_MMSAWT"], g["_STSTR"], np.arange(len(g)))
        lo, hi = logit_ci(r["p"], r["se"])
        name = g["MMSANAME"].iloc[0]
        rows.append({"mmsa": str(int(code)), "mmsa_name": name.decode() if isinstance(name, bytes) else str(name),
                     "value": 100 * r["p"], "se": 100 * r["se"], "ci_low": 100 * lo, "ci_high": 100 * hi,
                     "n_unweighted": r["n"], "n_cases": r["cases"]})
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- audit
def _pct(n, d):
    return f"{n:,} ({100 * n / d:.2f}%)" if d else f"{n:,}"


def missingness(raw: pd.DataFrame, year: int) -> list[str]:
    lines = []
    n = len(raw)
    for v in DESIGN_VARS + COVARIATES + LC_VARS[year]:
        if v in raw:
            lines.append(f"| {year} | {v} | {_pct(int(raw[v].isna().sum()), n)} |")
    return lines


def write_audit(s: dict) -> Path:
    yr_lines = "\n".join(
        f"| {y} | {v['n_records']:,} | {v['n_states_50dc']} (+{v['n_territories']} territories) | "
        f"{', '.join(v['lc_vars_in_layout']) or 'none'} | {v['n_layout_vars']} | {v['substate_vars'] or 'none'} | "
        f"[codebook]({v['codebook_url']}) `{v['codebook_file']}` |" for y, v in s["years"].items())
    var_lines = "\n".join(f"| {y} | {var} | {e['section']} | {e['label']} | {e['question'][:220]} |"
                          for y, ents in s["codebook_lc"].items() for var, e in ents.items())
    cov_lines = "\n".join(f"| {y} | {var} | {e['label']} |" for y, ents in s["codebook_cov"].items()
                          for var, e in ents.items())
    a = s["agreement"]
    worst = "; ".join(f"{w['jurisdiction']} {w['ours_current_std']:.2f} vs {w['mmwr_current_pct']:.1f}"
                      for w in a["worst"])
    files = "\n".join(f"| {k} | {v['url'][:140]} | {v['bytes']:,} | `{v['sha256'][:16]}...` | {v['retrieved_at']} |"
                      for k, v in sorted(s["manifest"].items()))
    miss = "\n".join(s["missing_lines"])
    nat = s["national"]
    mm = s["mmsa"]
    text = f"""# DATA AUDIT — CDC BRFSS LLCP public data (Long COVID items, 2022-2024)

| Field | Value |
|---|---|
| source_id | {SOURCE_ID} |
| Source (dataset/API name, exact files/endpoints) | BRFSS annual combined landline + cell ("LLCP") public SAS transport files 2022, 2023, 2024; their codebooks (USCODE HTML in `codebook<yy>_llcp-v2-508.zip`) and one-column variable layouts; SMART BRFSS MMSA 2023 file; MMWR 2024;73(50) supplementary table (validation); NCHS 2013 urban-rural codes and OMB 2020/2023 CBSA delineations (county urbanicity and MMSA membership, used by `geography.sae`); ACS 2020-2024 Summary File tables {', '.join(ACS_POSTSTRAT_TABLES)} (post-stratification cells, used by `geography.sae`) |
| Publishing organization | CDC, Division of Population Health (BRFSS); U.S. Census Bureau (ACS, delineations); NCHS (urban-rural codes) |
| Retrieval date (UTC) | {s['retrieved_at']} |
| Source version / release | LLCP 2022 (codebook {s['years'][2022]['codebook_file']}), 2023 ({s['years'][2023]['codebook_file']}), 2024 ({s['years'][2024]['codebook_file']}) |
| Source update date / cadence | annual (each year's file released the following August/September) |
| License / access conditions | U.S. federal government work, public domain; open download, no registration |
| Unit of observation | adult respondent (person-level survey record); kept in `data/interim/cdc_brfss` only |
| Sample size (actual, as ingested) | see per-year table; 2023 analytic file {s['years'][2023]['n_records']:,} records, of which {s['n_2023_50dc']:,} in the 50 states + DC |
| Geography (resolution, vintage) | state (`_STATE`); county-derived urbanicity flags `_METSTAT` / `_URBSTAT`. **No county identifier on the public LLCP file** (layout scanned for county/ZIP/tract/MMSA variable names: none). MMSA (metropolitan area) only in the separate SMART MMSA file. |
| Person-level? | yes (survey microdata; interim only, no processed person table) |
| Geographic? | yes (state direct estimates) |
| Omics? | no |
| Wearable? | no |
| True participant linkage across modalities? | no — single cross-sectional survey; the LLCP and MMSA files are NOT linked record-to-record here |

## Years and Long COVID items (verified from layouts and codebooks)
| year | records | states (50+DC) | COVID variables in layout | layout variables | sub-state geography variables | codebook |
|---|---|---|---|---|---|---|
{yr_lines}

Codebook entries for the Long COVID items (question text verbatim from the codebook):

| year | variable | section | label | question |
|---|---|---|---|---|
{var_lines}

**2024: no Long COVID item** — the 2024 layout and codebook contain no `COVID*` variable, so 2024 is downloaded and
checked but contributes no Long COVID estimate. **2022 asks EVER** (COVIDSMP "Did you have any symptoms lasting 3
months or longer..."), **2023 asks CURRENT** (COVIDSM1 "Do you currently have symptoms..."); the two years are not
pooled. Current Long COVID (2023) is the measure modelled by `geography.sae`. Landing pages: {', '.join(LANDING.format(y=y) for y in YEARS)}.

Covariate and design variables (codebook labels):

| year | variable | label |
|---|---|---|
{cov_lines}

## Outcome definitions
* `lc_current_pct_all_adults` (primary for the SAE): 1 = COVIDPO1=1 and COVIDSM1=1; 0 = COVIDPO1=2, or COVIDPO1=1 and
  COVIDSM1=2; missing otherwise (don't know, refused, or blank). Blank COVIDPO1 occurs only in partial interviews
  (DISPCODE 1200) that ended before the COVID section ({s['blank_partial']}).
* `lc_current_pct_all_adults_mmwr_def`: the MMWR denominator (all adults; non-response counted as "no"). It reproduces the
  published numbers (below) but counts break-offs as non-cases; state break-off rates differ ({s['blank_range']}), so it is a
  state-varying downward bias. National crude: complete case {nat['cc']:.2f}% vs MMWR definition {nat['mmwr']:.2f}%.
* `lc_significant_limitation_pct_all_adults`: current Long COVID and COVIDACT=1 ("yes, a lot"), complete case.
* `lc_significant_limitation_pct_current_lc_mmwr_def`: COVIDACT=1 among COVIDSM1=1 (MMWR second column).
* `lc_ever_pct_all_adults` (2022): COVIDPOS in (1, 3) and COVIDSMP=1; complete case.

## Validation against the published MMWR jurisdiction table (2023)
Source: Ford ND et al., "Notes from the Field: Long COVID and Significant Long COVID-Associated Activity Limitation Among
Adults, by Jurisdiction — United States, 2023", MMWR 2024;73(50) ({MMWR_URL}); supplementary table {MMWR_SUPP_URL}
(parsed with pdftotext; {a['n_jurisdictions']} jurisdictions with estimates). MMWR estimates are weighted, SUDAAN, and
standardised by sex x age (18-44, 45-64, 65+) to the 2020 Census civilian noninstitutionalised adults; ours are
standardised to ACS 2020-2024 adults (the exact MMWR standard population is not published in the table), Taylor-linearised
SEs (strata `_STSTR`, PSU `_PSU`), logit CIs.

* Unweighted case counts identical to MMWR in **{a['case_counts_identical']} / {a['n_jurisdictions']}** jurisdictions.
* Current Long COVID, standardised: mean absolute difference **{a['current_mae_pp']:.3f} percentage points**, max {a['current_max_abs_pp']:.2f} pp;
  {a['current_within_0_1pp']} / {a['n_jurisdictions']} within 0.1 pp (the table's rounding), {a['current_within_0_3pp']} within 0.3 pp; Pearson r = {a['current_pearson']:.4f}.
  Largest differences: {worst}. Crude (unstandardised) MAE {a['current_mae_crude_pp']:.3f} pp.
* National: the MMWR 'National' row counts {a['national_cases_mmwr']:,} cases (all jurisdictions incl. Guam, Puerto Rico,
  US Virgin Islands); ours is the 50 states + DC ({a['national_cases_ours']:,} cases): standardised {a['national_ours_std']:.2f}%, crude {a['national_ours_crude']:.2f}% vs MMWR {a['national_mmwr']:.1f}%.
* Significant limitation among adults with current Long COVID: MAE {a['sig_mae_pp']:.3f} pp, max {a['sig_max_abs_pp']:.2f} pp, r = {a['sig_pearson']:.4f}.
* Full comparison: `results/tables/brfss_mmwr_2023_comparison.csv`.

## SMART BRFSS MMSA file (2023)
The MMSA file carries COVIDPO1 / COVIDSM1 / COVIDACT, `_MMSA` and `_MMSAWT` for {mm['n_mmsa']} metropolitan/micropolitan areas
or divisions (>= 500 respondents; OMB March 2020 delineation per the SMART methodology PDF), {mm['n_records']:,} records.
It has no `_PSU`; MMSA direct estimates treat respondents as PSUs within `_STSTR`. Used ONLY as an external sub-state
validation of the county SAE (`results/BRFSS_SAE_RESULTS.md`); the MMSA file is not linked to LLCP records, and MMSA
estimates are not written as burden rows. Table: `results/tables/brfss_mmsa_2023_direct.csv`.

## Files / endpoints retrieved
| file | url | bytes | sha256 | retrieved_at |
|---|---|---|---|---|
{files}

## Missingness (raw values, per year, all records)
| year | variable | missing |
|---|---|---|
{miss}

`_AGEG5YR` = 14 (don't know/refused) and `_EDUCAG` = 9 are treated as missing for modelling; `_IMPRACE` (imputed) has
no missing values by construction; `_INCOMG1` = 9 (unknown) affects ~20% and income is therefore not a model covariate.

## Linkage strategy
`_STATE` -> 2-character state FIPS. No county or ZIP is on the public file, so no respondent is placed in a county.
County estimates exist only as the model output of `measure_it.geography.sae` (post-stratification), documented with
its own validation. The urbanicity flags are respondent attributes derived by CDC from the respondent's county and
enter the model as individual covariates.

## Limitations and caveats
* Telephone survey, self-reported Long COVID, no clinical confirmation; median response rate 2023 44.7% (21.7-63.1%).
* Kentucky and Pennsylvania are absent in 2023 (CDC: insufficient data); their state values are model predictions only.
* 2022 (ever) and 2023 (current) measure different things; 2024 has no item.
* The MMWR definition counts partial-interview break-offs as "no"; the SAE primary uses complete cases.
* Weighted with `_LLCPWT` (raked to age, sex, race, education, marital status, tenure, phone ownership and, where
  available, region/county margins).

## Processed outputs
| table | rows |
|---|---|
| brfss_long_covid_direct_estimates | {s['n_direct']:,} |

## Reproduce
`uv run python -m measure_it.ingestion.brfss` (then `uv run python -m measure_it.geography.sae`)
"""
    p = raw_dir(SOURCE_ID) / "DATA_AUDIT.md"
    # keep the section written by measure_it.geography.sae (it documents the model built on this source)
    start, end = "<!-- SAE SECTION START", "<!-- SAE SECTION END -->"
    if p.exists() and start in (old := p.read_text()):
        text = text.rstrip("\n") + "\n\n" + start + old.split(start, 1)[1].split(end, 1)[0] + end + "\n"
    p.write_text(text)
    return p


# --------------------------------------------------------------------------- main
def load_year(y: int) -> tuple[pd.DataFrame, dict]:
    urls = year_urls(y)
    lay = download_file(urls["layout"][0], SOURCE_ID, urls["layout"][1], allow_html=True).read_text(errors="replace")
    layout_vars = layout_variables(lay)
    cb_file, cb = read_codebook(download_file(urls["codebook"][0], SOURCE_ID, urls["codebook"][1]))
    lc_in_layout = sorted(v for v in layout_vars if v.startswith("COVID"))
    expected = LC_VARS[y]
    missing = [v for v in expected + DESIGN_VARS + COVARIATES if v not in layout_vars or v not in cb]
    if missing:
        raise ValueError(f"BRFSS {y}: expected variables not in layout/codebook: {missing}")
    lc_in_codebook = sorted(v for v in cb if v.startswith("COVID"))
    if not expected and (lc_in_layout or lc_in_codebook):
        raise ValueError(f"BRFSS {y} now lists COVID variables {lc_in_layout or lc_in_codebook}; review LC_VARS")
    raw = read_subset(download_file(urls["xpt"][0], SOURCE_ID, urls["xpt"][1], timeout=1800, max_retries=3),
                      DESIGN_VARS + COVARIATES + expected, interim_dir(SOURCE_ID) / f"llcp{y}_subset.parquet")
    info = {"n_records": len(raw), "lc_vars_in_layout": lc_in_layout, "n_layout_vars": len(layout_vars),
            "substate_vars": ", ".join(v for v in layout_vars if SUBSTATE_PATTERN.search(v)),
            "codebook_url": urls["codebook"][0], "codebook_file": cb_file,
            "codebook_lc": {v: cb[v] for v in expected},
            "codebook_cov": {v: cb[v] for v in DESIGN_VARS + COVARIATES}}
    return raw, info


def analysis_frame(year: int = 2023) -> pd.DataFrame:
    """Derived 2023 (or 2022) analysis frame from the interim subset (for geography.sae)."""
    raw, _ = load_year(year)
    return derive(raw, year)


def run() -> dict:
    fetch()
    names = state_names()
    universe = {f for f in names if int(f) <= 56}
    years, raws, derived = {}, {}, {}
    for y in YEARS:
        raw, info = load_year(y)
        st = raw["_STATE"].astype(int).astype(str).str.zfill(2)
        info["n_states_50dc"] = int(st[st.isin(universe)].nunique())
        info["n_territories"] = int(st[~st.isin(universe)].nunique())
        years[y], raws[y] = info, raw
        if LC_VARS[y]:
            derived[y] = derive(raw, y)
    standard = acs_national_standard()
    d23 = derived[2023]
    d23_us = d23[d23["state_fips"].isin(universe)]
    m23 = [m for m, v in MEASURES.items() if v[0] == 2023]
    m22 = [m for m, v in MEASURES.items() if v[0] == 2022]
    direct = pd.concat([direct_estimates(d23, m23, standard).assign(year=2023),
                        direct_estimates(derived[2022], m22, standard).assign(year=2022)], ignore_index=True)
    direct["in_50_states_dc"] = direct["geo_id"].isin(universe | {"US50DC"})

    mmwr_text = pdf_text(raw_dir(SOURCE_ID) / MMWR_SUPP_FILE)
    if mmwr_text is None:
        raise RuntimeError("pdftotext is required to parse the MMWR supplementary table")
    mmwr = parse_mmwr_table(mmwr_text)
    cmp_ = mmwr_comparison(direct[direct.year == 2023], mmwr)
    agreement = agreement_summary(cmp_)
    from ..config import TABLES
    TABLES.mkdir(parents=True, exist_ok=True)
    cmp_.to_csv(TABLES / "brfss_mmwr_2023_comparison.csv", index=False)
    mm = mmsa_direct(2023)
    mm.to_csv(TABLES / "brfss_mmsa_2023_direct.csv", index=False)

    manifest = load_manifest(SOURCE_ID)["files"]
    retrieved = manifest["2023/LLCP2023XPT.zip"]["retrieved_at"]
    version = "BRFSS LLCP 2022, 2023, 2024 public files (" + "; ".join(
        f"{y}: {years[y]['codebook_file']}" for y in YEARS) + ")"
    out = direct.copy()
    out["geo_id"] = out["geo_id"].replace({"US50DC": "US"})
    out["condition_id"] = "long_covid"
    out["measure_label"] = out["measure_id"].map(lambda m: MEASURES[m][1])
    out["denominator_population"] = out["measure_id"].map(lambda m: MEASURES[m][2])
    out["metric_type"] = "prevalence_pct"
    out["value_unit"] = "percent"
    out["ci_level"] = "95% (logit, Taylor-linearised)"
    out["burden_evidence_level"] = "A"
    out["primary_burden_measure"] = (out["measure_id"].eq("lc_current_pct_all_adults")
                                     & out["estimate_type"].eq("crude") & out["in_50_states_dc"])
    out["period_start"] = out["year"].astype(str) + "-01-01"
    out["period_end"] = out["year"].astype(str) + "-12-31"
    out = add_provenance(
        out, data_layer="geographic", source_name=SOURCE_NAME, source_version=version, retrieved_at=retrieved,
        evidence_type="survey_estimate",
        source_record_id=lambda d: ("brfss" + d["year"].astype(str) + ":" + d["measure_id"] + ":" + d["geo_id"]
                                    + ":" + d["estimate_type"]),
        source_geographic_resolution="state", evidence_level="A",
        provenance_notes="Direct design-based BRFSS estimate (weights _LLCPWT, strata _STSTR, PSU _PSU). State or "
                         "national only; never county prevalence.")
    out.loc[out["geo_level"] == "national", "source_geographic_resolution"] = "national"
    write_table(out, "brfss_long_covid_direct_estimates", producer=PRODUCER,
                description="BRFSS 2022 (ever) / 2023 (current, significant limitation) Long COVID direct estimates, "
                            "national + state")

    blank = raws[2023]["COVIDPO1"].isna()
    partial = raws[2023]["DISPCODE"].eq(1200)
    st23 = raws[2023]["_STATE"].astype(int).astype(str).str.zfill(2)
    rate = blank[st23.isin(universe)].groupby(st23[st23.isin(universe)]).mean()
    nat_rows = direct[(direct.geo_id == "US50DC") & (direct.estimate_type == "crude")].set_index("measure_id")["value"]
    stats = {
        "years": years, "codebook_lc": {y: years[y]["codebook_lc"] for y in YEARS if LC_VARS[y]},
        "codebook_cov": {2023: years[2023]["codebook_cov"]},
        "agreement": agreement, "manifest": manifest, "retrieved_at": retrieved,
        "missing_lines": sum((missingness(raws[y], y) for y in YEARS), []),
        "n_2023_50dc": int(len(d23_us)), "n_direct": len(out),
        "blank_partial": f"{int((blank & partial).sum()):,} of {int(blank.sum()):,} blank COVIDPO1 are partial interviews",
        "blank_range": f"blank COVIDPO1 {100 * rate.min():.1f}% to {100 * rate.max():.1f}% of records by state",
        "national": {"cc": nat_rows["lc_current_pct_all_adults"], "mmwr": nat_rows["lc_current_pct_all_adults_mmwr_def"]},
        "mmsa": {"n_mmsa": int(len(mm)), "n_records": int(mm["n_unweighted"].sum())},
    }
    write_audit(stats)
    write_registry_entry({
        "source_id": SOURCE_ID, "name": "CDC Behavioral Risk Factor Surveillance System (BRFSS) LLCP public data, "
                                        "Long COVID items 2022-2023 (2024 checked: no item)",
        "publisher": "CDC Division of Population Health", "landing_url": LANDING.format(y=2023),
        "access_urls": [year_urls(y)["xpt"][0] for y in YEARS] + [year_urls(y)["codebook"][0] for y in YEARS]
                       + [year_urls(2023)["mmsa_xpt"][0], MMWR_SUPP_URL],
        "license": "U.S. federal government work (public domain)", "access_conditions": "open download, no registration",
        "retrieved_at": retrieved, "source_version": version, "update_date": "annual",
        "data_layer": "geographic",
        "unit_of_observation": "adult survey respondent (interim only); processed = state/national direct estimates",
        "sample_size": {"records_2022": years[2022]["n_records"], "records_2023": years[2023]["n_records"],
                        "records_2024": years[2024]["n_records"], "records_2023_50_states_dc": int(len(d23_us)),
                        "current_lc_cases_2023": int(d23["lc_current_pct_all_adults_mmwr_def"].sum()),
                        "mmsa_2023_areas": int(len(mm))},
        "geographic_resolution": "state (public LLCP); urbanicity flags; MMSA in separate SMART file; no county",
        "person_level": True, "geographic": True, "omics": False, "wearable": False,
        "participant_linkage": "none (cross-sectional survey; LLCP and MMSA files not linked)",
        "true_participant_linkage_across_modalities": False, "status": "ingested",
        "processed_outputs": ["brfss_long_covid_direct_estimates", "geo_condition_burden__brfss_long_covid_sae",
                              "acs_county_adult_poststrat_cells"],
        "audit": f"data/raw/{SOURCE_ID}/DATA_AUDIT.md", "ingestion_module": PRODUCER,
        "burden_evidence_level": "A", "condition_ids": ["long_covid"],
        "validation": {"mmwr_2023_state_mae_pp": round(agreement["current_mae_pp"], 4),
                       "mmwr_2023_case_counts_identical": f"{agreement['case_counts_identical']}/{agreement['n_jurisdictions']}"},
        "limitations": [
            "Self-reported Long COVID; telephone survey with ~45% median response rate",
            "No county identifier on the public file; county values only via geography.sae (model)",
            "Kentucky and Pennsylvania absent from 2023 public data",
            "2022 item is EVER long COVID, 2023 is CURRENT; 2024 has no Long COVID item",
            "MMWR definition counts break-offs as 'no' (downward, state-varying bias); SAE uses complete cases",
        ],
    })
    summary = {"records": {y: years[y]["n_records"] for y in YEARS}, "agreement": {k: v for k, v in agreement.items()
                                                                                   if k != "worst"},
               "national": stats["national"], "mmsa": stats["mmsa"]}
    print(json.dumps(summary, indent=2, default=str))
    return summary


if __name__ == "__main__":
    run()
