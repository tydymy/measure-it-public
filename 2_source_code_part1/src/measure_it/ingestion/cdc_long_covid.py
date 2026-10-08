"""CDC Long COVID prevalence: NCHS "Post-COVID Conditions" (Household Pulse Survey), national + state.

Source: data.cdc.gov dataset gsea-w83j ("Post-COVID Conditions", NCHS; Household Pulse Survey
Phases 3.5-4.2, 2022-06-01 .. 2024-09-16). Values are survey-weighted percentages of adults 18+
with 95% CIs. States have one row per indicator x period; demographic subgroups are national only.

Successor check (recorded on every run): the Census Household Trends and Outlook Pulse Survey
(HTOPS), into which the HPS moved in 2025, is scanned for long-COVID items by downloading its
published questionnaires and searching their text (pdftotext). BRFSS 2022-2023 carried core
long-COVID questions (COVIDSMP / COVIDSM1 / COVIDACT), but CDC's BRFSS prevalence datasets on
data.cdc.gov (dttw-5yxu, d2rk-yvas) publish no COVID indicator; the check queries them. No
state estimate is derived here from BRFSS microdata.

Guardrail: geography stays at the source resolution (national, state). Never county. A state
value copied onto counties by a downstream step must keep source_geographic_resolution='state'.

Output (data/processed): geo_condition_burden__cdc_long_covid
  condition_id=long_covid, burden_evidence_level='A' (direct self-reported condition measure).
  Indicators are kept separate (ever / current / activity limitation, each with its denominator).
  'Ever had COVID' (infection history) is not a long-COVID burden measure and is not written.

Run: uv run python -m measure_it.ingestion.cdc_long_covid
"""
from __future__ import annotations

import json
import re
import shutil
import subprocess
from pathlib import Path

import pandas as pd

from ..config import raw_dir
from ..download import download_file, load_manifest
from ..geography.crosswalk import STATE_FIPS
from ..http import get_json
from ..provenance import add_provenance
from ..registry import write_registry_entry
from ..store import read_table, write_table

SOURCE_ID = "cdc_long_covid"
SOURCE_NAME = "CDC NCHS Post-COVID Conditions (Household Pulse Survey)"
PRODUCER = "measure_it.ingestion.cdc_long_covid"
DOMAIN = "data.cdc.gov"
DATASET_ID = "gsea-w83j"
LANDING_URL = "https://www.cdc.gov/nchs/covid19/pulse/long-covid.htm"
HTOPS_QUESTIONNAIRE_PAGES = [
    "https://www.census.gov/programs-surveys/household-pulse-survey/technical-documentation/questionnaires.2025.html",
    "https://www.census.gov/programs-surveys/household-pulse-survey/technical-documentation/questionnaires.2026.html",
]
BRFSS_PREVALENCE_DATASETS = {"dttw-5yxu": "BRFSS Prevalence Data (2011 to present)",
                             "d2rk-yvas": "BRFSS Age-Adjusted Prevalence Data (2011 to present)"}
BRFSS_QUESTIONS_DATASET = "iuq5-y9ct"
NHIS_DATASETS = ["25m4-6qqq", "krhz-spsc", "trpk-sp8z", "wpti-gvdi"]
LONG_COVID_TEXT = re.compile(r"long[- ]covid|post[- ]covid|symptoms lasting (?:3|three) months|"
                             r"(?:3|three) months or longer", re.IGNORECASE)

# verbatim indicator -> (measure_id, indicator_family, denominator population)
INDICATORS: dict[str, tuple[str, str, str]] = {
    "Ever experienced long COVID, as a percentage of all adults":
        ("lc_ever_pct_all_adults", "ever_had_long_covid", "all adults 18+"),
    "Ever experienced long COVID, as a percentage of adults who ever had COVID":
        ("lc_ever_pct_adults_ever_covid", "ever_had_long_covid", "adults 18+ who ever had COVID"),
    "Currently experiencing long COVID, as a percentage of all adults":
        ("lc_current_pct_all_adults", "currently_has_long_covid", "all adults 18+"),
    "Currently experiencing long COVID, as a percentage of adults who ever had COVID":
        ("lc_current_pct_adults_ever_covid", "currently_has_long_covid", "adults 18+ who ever had COVID"),
    "Any activity limitations from long COVID, as a percentage of all adults":
        ("lc_any_activity_limitation_pct_all_adults", "activity_limitation_any", "all adults 18+"),
    "Any activity limitations from long COVID, as a percentage of adults who currently have long COVID":
        ("lc_any_activity_limitation_pct_current_lc", "activity_limitation_any", "adults 18+ currently with long COVID"),
    "Significant activity limitations from long COVID, as a percentage of all adults":
        ("lc_significant_activity_limitation_pct_all_adults", "activity_limitation_significant", "all adults 18+"),
    "Significant activity limitations from long COVID, as a percentage of adults who currently have long COVID":
        ("lc_significant_activity_limitation_pct_current_lc", "activity_limitation_significant",
         "adults 18+ currently with long COVID"),
}
EXCLUDED_INDICATORS = {"Ever had COVID": "COVID infection history (context/denominator), not a long-COVID burden measure"}
PRIMARY_MEASURE = "lc_current_pct_all_adults"


# --------------------------------------------------------------------------- pure helpers
def indicator_map(indicator: str) -> tuple[str, str, str] | None:
    """Verbatim indicator -> (measure_id, family, denominator); None for excluded; raises if unknown."""
    if indicator in INDICATORS:
        return INDICATORS[indicator]
    if indicator in EXCLUDED_INDICATORS:
        return None
    raise ValueError(f"unknown Post-COVID Conditions indicator {indicator!r}; review INDICATORS before ingesting")


def hps_phase_label(phase: str, start_mmddyyyy: str) -> str:
    """The raw Phase column is numeric, so HPS Phase 3.10 (Aug-Oct 2023) is stored as 3.1.

    Long-COVID items start in Phase 3.5 (June 2022); HPS Phase 3.1 (2021) predates them, so any
    '3.1' dated after 2022 is Phase 3.10.
    """
    if str(phase) == "3.1" and pd.to_datetime(start_mmddyyyy, format="%m/%d/%Y") >= pd.Timestamp("2022-06-01"):
        return "3.10"
    return str(phase)


def state_name_to_fips(geos: pd.DataFrame) -> dict[str, str]:
    st = geos[geos["geo_level"] == "state"]
    return dict(zip(st["name"], st["geo_id"]))


def tidy(raw: pd.DataFrame, name_to_fips: dict[str, str]) -> tuple[pd.DataFrame, dict]:
    """Raw gsea-w83j CSV -> long burden table. Returns (table, counts of dropped rows)."""
    d = raw.copy()
    dropped = {"collection_gap_rows": int((d["Phase"] == "-1").sum())}
    d = d[d["Phase"] != "-1"]  # placeholder rows for weeks between HPS phases: no data were collected
    mapped = d["Indicator"].map(indicator_map)
    dropped["ever_had_covid_rows"] = int(mapped.isna().sum())
    d, mapped = d[mapped.notna()], mapped[mapped.notna()]
    national = d["State"].eq("United States")
    out = pd.DataFrame({
        "geo_id": d["State"].map(name_to_fips).where(~national, "US"),
        "geo_level": national.map({True: "national", False: "state"}),
        "geo_name": d["State"],
        "condition_id": "long_covid",
        "measure_id": [m[0] for m in mapped],
        "measure_label": d["Indicator"],
        "indicator_family": [m[1] for m in mapped],
        "denominator_population": [m[2] for m in mapped],
        "metric_type": "prevalence_pct",
        "value": pd.to_numeric(d["Value"], errors="coerce"),
        "value_unit": "percent",
        "ci_low": pd.to_numeric(d["LowCI"], errors="coerce"),
        "ci_high": pd.to_numeric(d["HighCI"], errors="coerce"),
        "ci_level": "95%",
        "subgroup_type": d["Group"],
        "subgroup": d["Subgroup"],
        "hps_phase_raw": d["Phase"],
        "hps_phase": [hps_phase_label(ph, st) for ph, st in zip(d["Phase"], d["Time Period Start Date"])],
        "time_period_id": pd.to_numeric(d["Time Period"], errors="coerce").astype("Int64"),
        "period_label": d["Time Period Label"],
        "period_start": pd.to_datetime(d["Time Period Start Date"], format="%m/%d/%Y").dt.date.astype(str),
        "period_end": pd.to_datetime(d["Time Period End Date"], format="%m/%d/%Y").dt.date.astype(str),
        "quartile_range": d["Quartile range"],
        "quartile_number": pd.to_numeric(d["Quartile number"], errors="coerce").astype("Int64"),
        "suppressed": d["Suppression Flag"].eq("1"),
    })
    out["year"] = pd.to_datetime(out["period_start"]).dt.year.astype("Int64")
    is_state = out["geo_level"] == "state"
    out["state_fips"] = out["geo_id"].where(is_state)
    out["state_abbr"] = out["state_fips"].map(STATE_FIPS)
    out["state_name"] = out["geo_name"].where(is_state)
    out["county_fips"] = pd.Series(pd.NA, index=out.index, dtype="str")
    out["geo_vintage"] = is_state.map({True: "2024", False: "national"})
    # national rows (geo_id 'US') have no row in `geographies`; same convention as geo_context__cdc_places
    out["in_canonical_geographies"] = out["geo_id"].isin(set(name_to_fips.values())) & is_state
    out["numerator"] = pd.array([pd.NA] * len(out), dtype="Int64")
    out["denominator"] = pd.array([pd.NA] * len(out), dtype="Int64")
    out["denominator_source"] = "survey-weighted percent; weighted denominators not published in gsea-w83j"
    out["primary_burden_measure"] = out["measure_id"].eq(PRIMARY_MEASURE) & out["subgroup_type"].isin(
        ["National Estimate", "By State"])
    out["burden_evidence_level"] = "A"
    out["survey"] = "Household Pulse Survey (Census Bureau + NCHS), experimental"
    if out.loc[out["geo_level"] == "state", "geo_id"].isna().any():
        bad = out.loc[(out["geo_level"] == "state") & out["geo_id"].isna(), "geo_name"].unique()
        raise ValueError(f"state names without FIPS: {list(bad)}")
    return out.reset_index(drop=True), dropped


def period_lengths(table: pd.DataFrame) -> str:
    p = table[["hps_phase", "period_start", "period_end"]].drop_duplicates()
    days = (pd.to_datetime(p["period_end"]) - pd.to_datetime(p["period_start"])).dt.days + 1
    major = p["hps_phase"].str.split(".").str[0]
    parts = [f"Phase {m}.x {int(days[major == m].min())}-{int(days[major == m].max())} days" for m in sorted(major.unique())]
    return "; ".join(parts)


def scan_text_for_long_covid(text: str) -> dict:
    hits = LONG_COVID_TEXT.findall(text)
    return {"long_covid_hits": len(hits), "covid_mentions": len(re.findall(r"covid", text, re.IGNORECASE)),
            "chars": len(text)}


def pdf_text(path: Path) -> str | None:
    exe = shutil.which("pdftotext")
    if not exe:
        return None
    r = subprocess.run([exe, "-q", str(path), "-"], capture_output=True, timeout=120, check=False)
    return r.stdout.decode("utf-8", errors="replace")


# --------------------------------------------------------------------------- successor checks
def check_htops() -> list[dict]:
    """Download every HTOPS questionnaire linked from the Census pages and scan for long-COVID items."""
    results = []
    qdir = raw_dir(SOURCE_ID) / "htops_questionnaires"
    qdir.mkdir(parents=True, exist_ok=True)
    for page in HTOPS_QUESTIONNAIRE_PAGES:
        try:
            html = download_file(page, SOURCE_ID, f"htops_questionnaires/{page.rsplit('/', 1)[-1]}",
                                 allow_html=True).read_text(errors="replace")
        except Exception as exc:  # noqa: BLE001 - record and continue
            results.append({"page": page, "error": str(exc)})
            continue
        for url in sorted(set(re.findall(r'href="(https://www2\.census\.gov/[^"]*HTOPS_[^"]*_ENGLISH\.pdf)"', html))):
            fname = url.rsplit("/", 1)[-1]
            try:
                p = download_file(url, SOURCE_ID, f"htops_questionnaires/{fname}")
            except Exception as exc:  # noqa: BLE001
                results.append({"page": page, "url": url, "error": str(exc)})
                continue
            text = pdf_text(p)
            row = {"page": page, "url": url, "file": fname}
            row.update(scan_text_for_long_covid(text) if text is not None else {"scan": "pdftotext not available"})
            results.append(row)
    return results


def check_brfss() -> dict:
    out = {}
    for ds, name in BRFSS_PREVALENCE_DATASETS.items():
        hits = get_json(f"https://{DOMAIN}/resource/{ds}.json", params={
            "$select": "question,count(*)", "$group": "question", "$limit": 50,
            "$where": "upper(question) like '%COVID%' OR upper(topic) like '%COVID%' OR upper(class) like '%COVID%'"})
        yr = get_json(f"https://{DOMAIN}/resource/{ds}.json", params={"$select": "min(year),max(year)"})[0]
        out[ds] = {"name": name, "covid_questions": len(hits), "years": f"{yr.get('min_year')}-{yr.get('max_year')}"}
    q = get_json(f"https://{DOMAIN}/resource/{BRFSS_QUESTIONS_DATASET}.json", params={
        "$where": "upper(topic) like '%COVID EFFECTS%'", "$limit": 50})
    out["brfss_questionnaire_items"] = sorted({f"{r['year']} {r['variablename']} ({r['type']})" for r in q})
    nhis = {}
    for ds in NHIS_DATASETS:
        rows = get_json(f"https://{DOMAIN}/resource/{ds}.json", params={
            "$select": "outcome_or_indicator,count(*)", "$group": "outcome_or_indicator", "$limit": 1000})
        names = [r["outcome_or_indicator"] for r in rows]
        nhis[ds] = {"indicators": len(names), "covid_indicators": [n for n in names if re.search(r"covid", n, re.IGNORECASE)]}
    out["nhis"] = nhis
    return out


# --------------------------------------------------------------------------- audit
def _pct(n: int, d: int) -> str:
    return f"{n:,} ({100 * n / d:.2f}%)" if d else f"{n:,}"


def write_audit(s: dict) -> Path:
    t = s["table"]
    ind = (t.groupby(["measure_id", "denominator_population"])
           .agg(rows=("value", "size"), suppressed=("suppressed", "sum"),
                first=("period_start", "min"), last=("period_end", "max"),
                periods=("time_period_id", "nunique")).reset_index())
    ind_md = "\n".join(f"| {r.measure_id} | {r.denominator_population} | {r.rows:,} | {_pct(int(r.suppressed), r.rows)} | "
                       f"{r.periods} | {r.first} .. {r.last} |" for r in ind.itertuples())
    grp = t.groupby("subgroup_type").agg(rows=("value", "size"), suppressed=("suppressed", "sum"),
                                         geos=("geo_id", "nunique"), subgroups=("subgroup", "nunique")).reset_index()
    grp_md = "\n".join(f"| {r.subgroup_type} | {r.rows:,} | {r.geos} | {r.subgroups} | {_pct(int(r.suppressed), r.rows)} |"
                       for r in grp.itertuples())
    miss = "\n".join(f"| {c} | {_pct(int(t[c].isna().sum()), len(t))} |" for c in
                     ["value", "ci_low", "ci_high", "quartile_number", "geo_id", "period_start"])
    ht = s["htops"]
    ht_md = "\n".join(
        f"| {r.get('file', r.get('page'))} | {r.get('long_covid_hits', r.get('scan', r.get('error')))} | "
        f"{r.get('covid_mentions', '')} |" for r in ht)
    br = s["brfss"]
    files = "\n".join(f"| {k} | {v['url'][:150]} | {v['bytes']:,} | `{v['sha256'][:16]}...` | {v['retrieved_at']} |"
                      for k, v in s["manifest"].items())
    nat = t[(t["measure_id"] == PRIMARY_MEASURE) & (t["subgroup_type"] == "National Estimate")].sort_values("period_start")
    first, last = nat.iloc[0], nat.iloc[-1]
    text = f"""# DATA AUDIT — CDC Long COVID prevalence (Household Pulse Survey "Post-COVID Conditions")

| Field | Value |
|---|---|
| source_id | cdc_long_covid |
| Source (dataset/API name, exact files/endpoints) | NCHS "Post-COVID Conditions", data.cdc.gov `{DATASET_ID}` (full CSV export). Successor checks: Census HTOPS questionnaires; BRFSS prevalence datasets `dttw-5yxu`, `d2rk-yvas`; BRFSS question list `{BRFSS_QUESTIONS_DATASET}` |
| Publishing organization | CDC National Center for Health Statistics with the U.S. Census Bureau (Household Pulse Survey) |
| Retrieval date (UTC) | {s['retrieved_at']} |
| Source version / release | Socrata rowsUpdatedAt {s['rows_updated']}; metadata "Temporal Applicability" {s['temporal']}; HPS phases {s['phases']} |
| Source update date / cadence | Final. Last data period ends {s['last_period']}; the HPS long-COVID series was not continued (see successor check) |
| License / access conditions | {s['license']}; open download |
| Unit of observation | indicator x geography (US or state) x subgroup x HPS collection period; survey-weighted % of adults 18+ with 95% CI |
| Sample size (actual, as ingested) | raw {s['raw_rows']:,} rows; written {len(t):,} rows ({s['dropped']['collection_gap_rows']} between-phase placeholder rows with no data and {s['dropped']['ever_had_covid_rows']:,} 'Ever had COVID' rows not written); {t['time_period_id'].nunique()} collection periods; {t.loc[t.geo_level == 'state', 'geo_id'].nunique()} states (incl. DC) + national. HPS respondent counts are not published in this dataset. |
| Geography (resolution, vintage) | national and state (50 states + DC). Subgroups (age, sex, race/ethnicity, education, gender identity, sexual orientation, disability) are national only. **No county data exist.** |
| Person-level? | no |
| Geographic? | yes |
| Omics? | no |
| Wearable? | no |
| True participant linkage across modalities? | no — published aggregate survey estimates |

## Burden labelling
`condition_id = long_covid`, `burden_evidence_level = A` (direct self-reported condition measure; dataset description:
"These adults had COVID and had some symptoms that lasted three months or longer"), `metric_type = prevalence_pct`,
`source_geographic_resolution` = `state` or `national` per row. `primary_burden_measure = True` marks
`{PRIMARY_MEASURE}` for the national estimate and each state. National, first vs last period:
{first['value']}% ({first['ci_low']}-{first['ci_high']}) in {first['period_label']} vs {last['value']}% ({last['ci_low']}-{last['ci_high']}) in {last['period_label']}.

## Indicators (kept separate)
| measure_id | denominator | rows | suppressed | periods | first .. last |
|---|---|---|---|---|---|
{ind_md}

Activity-limitation indicators start {s['act_start']} (HPS phase {s['act_phase']}); the ever/current indicators start {s['lc_start']}.

## Groups
| subgroup_type | rows | geographies | subgroups | suppressed |
|---|---|---|---|---|
{grp_md}

## Is there anything after September 2024?
* **Census HTOPS (the 2025- successor of the Household Pulse Survey):** every English questionnaire linked from the
  2025 and 2026 questionnaire pages was downloaded and its text scanned (pdftotext) for
  `{LONG_COVID_TEXT.pattern}`:

| questionnaire | long-COVID pattern hits | 'COVID' mentions (vaccine item only) |
|---|---|---|
{ht_md}

  Result: **{s['htops_answer']}** The HPS long-COVID series ends with the {s['last_period']} period.
* **BRFSS:** BRFSS carried core long-COVID questions ({'; '.join(br['brfss_questionnaire_items'])}; source `{BRFSS_QUESTIONS_DATASET}`),
  but CDC's BRFSS prevalence datasets publish no COVID indicator: `dttw-5yxu` ({br['dttw-5yxu']['years']}) {br['dttw-5yxu']['covid_questions']} COVID questions;
  `d2rk-yvas` ({br['d2rk-yvas']['years']}) {br['d2rk-yvas']['covid_questions']}. State long-COVID prevalence from BRFSS therefore exists only in
  journal/MMWR tables or would have to be computed from BRFSS microdata with its survey design; neither is done here.
* **NHIS:** NHIS summary-statistics datasets on data.cdc.gov, long-COVID/COVID indicators found: {s['nhis_md']}.
  The NCHS long-COVID page states that long-COVID questions "were also included on the National Health Interview
  Survey (NHIS) in 2022"; no NHIS long-COVID
  estimate is published as a dataset among those checked.

## Files / endpoints retrieved
| file | url | bytes | sha256 | retrieved_at |
|---|---|---|---|---|
{files}

## Key variables
| column | raw | meaning |
|---|---|---|
| measure_id / measure_label | Indicator | one of 8 long-COVID indicators (ever / current / any or significant activity limitation), each with its own denominator |
| denominator_population | derived from Indicator | all adults, adults who ever had COVID, or adults currently with long COVID |
| value, ci_low, ci_high | Value, LowCI, HighCI | weighted percent, 95% CI |
| subgroup_type, subgroup | Group, Subgroup | National Estimate / By State / national demographic subgroups |
| time_period_id, period_label, period_start, period_end, hps_phase | Time Period*, Phase | HPS collection window (2-4 weeks) |
| quartile_range, quartile_number | Quartile* | CDC state map quartiles (state rows) |
| suppressed | Suppression Flag | 1 = estimate suppressed by NCHS (value blank) |

## Missingness
{_pct(int(t['suppressed'].sum()), len(t))} rows are suppressed (value, CI blank); every blank value is a suppressed row.

| column | missing |
|---|---|
{miss}

## Linkage strategy
`geo_id` = 2-char state FIPS (from `geographies` state names) or `US`. {s['state_rows_match']:,} of {s['state_rows']:,} state rows
match a `geographies` state row (Connecticut = state 09, {s['ct_rows']} rows; no county vintage question arises). Joins to county tables are allowed only
as inherited state context with `source_geographic_resolution='state'` kept; never as county prevalence.
Subgroup rows are national and must not be combined with state rows as if they were state subgroups.

## Limitations and caveats
* Household Pulse Survey is an experimental, online, cross-sectional survey with low response rates
  (Census nonresponse-bias report linked from the NCHS page); estimates may carry nonresponse bias.
* Self-reported long COVID with no clinical confirmation; definition = symptoms lasting 3+ months after COVID.
* Collection periods differ in length: {s['period_lengths']}.
* State estimates have wide CIs; many state x indicator x period cells are suppressed.
* Series ends {s['last_period']}; there is no public post-2024 state long-COVID series from CDC as data.

## Processed outputs
| table | rows |
|---|---|
| geo_condition_burden__cdc_long_covid | {len(t):,} |

## Reproduce
`uv run python -m measure_it.ingestion.cdc_long_covid`
"""
    p = raw_dir(SOURCE_ID) / "DATA_AUDIT.md"
    p.write_text(text)
    return p


# --------------------------------------------------------------------------- main
def run() -> dict:
    meta_path = download_file(f"https://{DOMAIN}/api/views/{DATASET_ID}.json", SOURCE_ID, f"metadata_{DATASET_ID}.json")
    meta = json.loads(meta_path.read_text())
    csv_path = download_file(f"https://{DOMAIN}/api/views/{DATASET_ID}/rows.csv?accessType=DOWNLOAD", SOURCE_ID,
                             f"post_covid_conditions_{DATASET_ID}.csv", min_bytes=100_000)
    raw = pd.read_csv(csv_path, dtype=str)
    geos = read_table("geographies")
    table, dropped = tidy(raw, state_name_to_fips(geos))

    manifest = load_manifest(SOURCE_ID)["files"]
    retrieved = manifest[csv_path.name]["retrieved_at"]
    cc = meta.get("metadata", {}).get("custom_fields", {}).get("Common Core", {})
    rows_updated = pd.Timestamp(meta["rowsUpdatedAt"], unit="s").isoformat() if meta.get("rowsUpdatedAt") else "UNKNOWN"
    version = f"data.cdc.gov {DATASET_ID}, rows updated {rows_updated}; HPS phases 3.5-4.2"
    out = add_provenance(
        table, data_layer="geographic", source_name=SOURCE_NAME, source_version=version, retrieved_at=retrieved,
        evidence_type="survey_estimate",
        source_record_id=lambda d: (f"{DATASET_ID}:" + d["measure_id"] + ":" + d["geo_name"] + ":"
                                    + d["subgroup_type"] + ":" + d["subgroup"] + ":" + d["time_period_id"].astype(str)),
        source_geographic_resolution="state", evidence_level="A",
        provenance_notes="Direct self-reported long-COVID measure (HPS). Burden level A at state/national resolution "
                         "only; never county prevalence. Subgroups are national.")
    out.loc[out["geo_level"] == "national", "source_geographic_resolution"] = "national"
    write_table(out, "geo_condition_burden__cdc_long_covid", producer=PRODUCER,
                description="HPS long-COVID prevalence (8 indicators), national + state, 2022-06..2024-09, level A")

    htops = check_htops()
    scanned = [r for r in htops if "long_covid_hits" in r]
    if scanned and all(r["long_covid_hits"] == 0 for r in scanned):
        htops_answer = f"no long-COVID item in any of the {len(scanned)} HTOPS questionnaires scanned."
    elif scanned:
        htops_answer = "long-COVID text FOUND in: " + ", ".join(r["file"] for r in scanned if r["long_covid_hits"])
    else:
        htops_answer = "not checked (questionnaires unavailable or pdftotext missing)."
    brfss = check_brfss()
    manifest = load_manifest(SOURCE_ID)["files"]
    stats = {
        "table": table, "dropped": dropped, "raw_rows": len(raw), "retrieved_at": retrieved,
        "rows_updated": rows_updated, "temporal": cc.get("Temporal Applicability", "UNKNOWN"),
        "phases": ", ".join(table.sort_values("period_start")["hps_phase"].unique()),
        "act_start": table.loc[table["indicator_family"].str.startswith("activity"), "period_start"].min(),
        "act_phase": table.loc[table["indicator_family"].str.startswith("activity")].sort_values("period_start")["hps_phase"].iloc[0],
        "lc_start": table.loc[~table["indicator_family"].str.startswith("activity"), "period_start"].min(),
        "period_lengths": period_lengths(table),
        "nhis_md": "; ".join(f"`{k}` {len(v['covid_indicators'])} of {v['indicators']}" for k, v in brfss["nhis"].items()),
        "last_period": table["period_end"].max(), "license": (meta.get("license") or {}).get("name", "UNKNOWN"),
        "htops": htops, "htops_answer": htops_answer, "brfss": brfss, "manifest": manifest,
        "state_rows": int((table["geo_level"] == "state").sum()),
        "state_rows_match": int(table.loc[table["geo_level"] == "state", "geo_id"].isin(
            geos.loc[geos["geo_level"] == "state", "geo_id"]).sum()),
        "ct_rows": int((table["geo_id"] == "09").sum()),
    }
    write_audit(stats)
    write_registry_entry({
        "source_id": SOURCE_ID, "name": "CDC NCHS Post-COVID Conditions (Household Pulse Survey long-COVID estimates)",
        "publisher": "CDC National Center for Health Statistics / U.S. Census Bureau", "landing_url": LANDING_URL,
        "access_urls": [f"https://{DOMAIN}/d/{DATASET_ID}", *HTOPS_QUESTIONNAIRE_PAGES,
                        *(f"https://{DOMAIN}/d/{k}" for k in BRFSS_PREVALENCE_DATASETS)],
        "license": stats["license"], "access_conditions": "open download, no registration",
        "retrieved_at": retrieved, "source_version": version,
        "update_date": f"rows updated {rows_updated}; series final (HPS long-COVID items ended {stats['last_period']})",
        "data_layer": "geographic",
        "unit_of_observation": "indicator x state/national x subgroup x HPS collection period (weighted % of adults)",
        "sample_size": {"rows": len(table), "raw_rows": len(raw), "periods": int(table["time_period_id"].nunique()),
                        "states_incl_dc": int(table.loc[table.geo_level == "state", "geo_id"].nunique()),
                        "indicators": int(table["measure_id"].nunique()),
                        "suppressed_rows": int(table["suppressed"].sum()),
                        "hps_respondents": "UNKNOWN / NOT AVAILABLE (not published in this dataset)"},
        "geographic_resolution": "national, state (no county)", "person_level": False, "geographic": True,
        "omics": False, "wearable": False, "participant_linkage": "not applicable (aggregate survey estimates)",
        "true_participant_linkage_across_modalities": False, "status": "ingested",
        "processed_outputs": ["geo_condition_burden__cdc_long_covid"],
        "audit": f"data/raw/{SOURCE_ID}/DATA_AUDIT.md", "ingestion_module": PRODUCER,
        "burden_evidence_level": "A", "condition_ids": ["long_covid"],
        "successor_check": {"htops": htops_answer,
                            "brfss_prevalence_covid_questions": {k: v["covid_questions"] for k, v in brfss.items()
                                                                 if k in BRFSS_PREVALENCE_DATASETS}},
        "limitations": [
            "Experimental online survey with low response rates; possible nonresponse bias",
            "Self-reported, no clinical confirmation",
            "State and national only; subgroups are national only; never county",
            "Series ends September 2024; HTOPS (2025-) has no long-COVID items",
            "Many state cells suppressed; wide state CIs",
        ],
    })
    summary = {"rows": len(table), "dropped": dropped, "htops": htops_answer,
               "states": int(table.loc[table.geo_level == "state", "geo_id"].nunique())}
    print(json.dumps(summary, indent=2, default=str))
    return summary


if __name__ == "__main__":
    run()
