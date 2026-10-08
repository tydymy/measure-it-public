"""`measure-it launch ...` commands (registered on the main CLI in measure_it.cli; docs/LAUNCH_AGENT.md).

    measure-it launch run <path> --condition long_covid --measurement nailfold_capillaroscopy [--out DIR]
                                 [--llm none|claude] [--approve] [--public] [--top 10] [--per-county 20] [--keep] ...
    measure-it launch profile <path> [--out DIR]
    measure-it launch map <path> --condition ... --measurement ... [--out DIR] [--llm none|claude]
    measure-it launch apply <run> [--approve]
    measure-it launch report <run>
"""
from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Optional

import typer

app = typer.Typer(no_args_is_help=True, add_completion=False,
                  help="Launch agent: any dataset -> harmonized layers (BYOD + All of Us-shaped OMOP) -> ranked, "
                       "traceable candidate launch sites (docs/LAUNCH_AGENT.md).")


def _out(msg: str) -> None:
    typer.echo(msg)
    sys.stdout.flush()


def _summary(res: dict) -> None:
    plan = res["plan"]
    _out(f"\nrun {plan['run']}: {plan['status']}" + (f"\n{plan['message']}" if plan.get("message") else ""))
    for s in plan["steps"]:
        why = s.get("reason") or s.get("error") or s.get("message") or ""
        _out(f"  {s['step']:<20} {s['status']:<8} {why[:140]}")
    la = res.get("launch")
    if la:
        prim = la["subgroup"] if la["primary"] == "subgroup" else la["engine"]
        _out(f"\n{la['definition']}")
        for c in (prim.get("counties") or [])[:10]:
            _out(f"  {c['rank']:>3}. {c['geo_name']} ({c['geo_id']})")
    if res.get("report"):
        _out(f"\nreport: {res['report'].get('report_md')}\nlaunch.json: {res['report'].get('launch_json')}")


def _data_use(deid: bool, permitted: bool, statement: Optional[str], ref: Optional[str]) -> dict | None:
    d = {"deidentified": True if deid else None, "use_permitted": True if permitted else None,
         "statement": statement, "irb_or_dua_reference": ref}
    return d if any(v is not None for v in d.values()) else None


@app.command("run")
def run_cmd(
    path: Path = typer.Argument(..., help="Input: a file or folder in (almost) any format, or a BYOD folder."),
    condition: str = typer.Option(..., "--condition", help="Condition id or set (e.g. long_covid, ptlds)."),
    measurement: str = typer.Option(..., "--measurement", help="Measurement bundle or class "
                                                               "(e.g. nailfold_capillaroscopy)."),
    out: Optional[Path] = typer.Option(None, "--out", help="Run folder (default launch_runs/<stamp>_<cond>_<meas>); "
                                                           "re-use it to resume after editing mapping.yaml."),
    llm: str = typer.Option("none", "--llm", help="none (default: no model is contacted) or claude (metadata-only "
                                                  "mapping proposals; needs ANTHROPIC_API_KEY)."),
    approve: bool = typer.Option(False, "--approve", help="Approve mapping.yaml (records approved_by / approved_at)."),
    public: bool = typer.Option(False, "--public", help="The input is a public, already de-identified dataset."),
    top: int = typer.Option(10, "--top", min=1, help="Top N counties."),
    per_county: int = typer.Option(20, "--per-county", min=1, help="Clinicians per county in the workbook."),
    keep: bool = typer.Option(False, "--keep", help="Leave the dataset ingested (default: byod remove at the end)."),
    synthetic: bool = typer.Option(False, "--synthetic", help="The input is generated test data (synthetic + demo)."),
    demo_rankings: bool = typer.Option(False, "--demo-rankings", help="Let a demo record enter the rankings for "
                                                                      "this run only (byod deploy --demo; restored "
                                                                      "by the teardown)."),
    confirm_deidentified: bool = typer.Option(False, "--confirm-deidentified",
                                              help="Confirm manifest data_use.deidentified: true."),
    confirm_use_permitted: bool = typer.Option(False, "--confirm-use-permitted",
                                               help="Confirm manifest data_use.use_permitted: true."),
    data_use_statement: Optional[str] = typer.Option(None, "--data-use-statement", help="data_use.statement."),
    irb_or_dua: Optional[str] = typer.Option(None, "--irb-or-dua", help="data_use.irb_or_dua_reference."),
    fast: bool = typer.Option(False, "--fast", help="Few CV repeats / permutations / bootstraps (smoke runs only)."),
    no_outreach: bool = typer.Option(False, "--no-outreach", help="Do not write the outreach workbook."),
    remap: bool = typer.Option(False, "--remap", help="Discard the run folder's mapping.yaml and propose again."),
    device_positive: Optional[str] = typer.Option(
        None, "--device-positive",
        help="Pre-specify the device-positive subgroup: a device column that is a 1/0 flag (e.g. cap_positive) or a "
             "threshold (e.g. 'cap_abnormality_score>=2'). Written into mapping.yaml before approval and locked by "
             "byod subgroup before any outcome is computed."),
    verbose: bool = typer.Option(False, "--verbose", help="Echo the engine log."),
) -> None:
    """Profile -> map -> approve -> apply (BYOD + OMOP) -> validate -> ingest -> evaluate -> subgroup -> similar ->
    deploy -> subgroup burden -> launch ranking -> outreach -> teardown -> report (plan.json records every step)."""
    from .orchestrate import run_launch
    try:
        res = run_launch(path, condition, measurement, out=out, llm=llm, approve=approve, public=public, top=top,
                         per_county=per_county, keep=keep, synthetic=synthetic, demo_rankings=demo_rankings,
                         data_use=_data_use(confirm_deidentified, confirm_use_permitted, data_use_statement,
                                            irb_or_dua),
                         fast=fast, outreach=not no_outreach, remap=remap, verbose=verbose, echo=_out,
                         device_positive=device_positive)
    except (ValueError, FileNotFoundError) as e:
        _out(f"launch refused: {e}")
        raise typer.Exit(2)
    _summary(res)
    st = res["plan"]["status"]
    raise typer.Exit(0 if st == "completed" else (3 if st == "stopped" else 1))


@app.command("profile")
def profile_cmd(path: Path = typer.Argument(...),
                out: Optional[Path] = typer.Option(None, "--out", help="Run folder for profile.json.")) -> None:
    """Profile the input (tables, columns, types, dictionary labels; metadata only) -> <out>/profile.json."""
    from .orchestrate import run_launch
    res = run_launch(path, "unspecified", "unspecified", out=out, outreach=False, stop_after="profile", echo=None)
    p = Path(res["out"]) / "profile.json"
    _summary(res)
    if p.exists():
        prof = json.loads(p.read_text())
        _out(f"\nprofile.json: {p}  ({prof.get('format')}, {len(prof.get('tables', []))} table(s))")
    raise typer.Exit(0 if p.exists() else 1)


@app.command("map")
def map_cmd(path: Path = typer.Argument(...),
            condition: str = typer.Option(..., "--condition"), measurement: str = typer.Option(..., "--measurement"),
            out: Optional[Path] = typer.Option(None, "--out"),
            llm: str = typer.Option("none", "--llm", help="none (default) or claude."),
            remap: bool = typer.Option(False, "--remap")) -> None:
    """Profile + mapping proposals -> <out>/mapping.yaml (human-editable; approve with `launch apply --approve`)."""
    from .orchestrate import run_launch
    res = run_launch(path, condition, measurement, out=out, llm=llm, outreach=False, remap=remap, stop_after="map",
                     echo=None)
    _summary(res)
    mp = Path(res["out"]) / "mapping.yaml"
    _out(f"\nmapping.yaml: {mp}\nnext: review it, then `measure-it launch apply {res['out']} --approve`")
    raise typer.Exit(0 if mp.exists() else 1)


@app.command("apply")
def apply_cmd(run: Path = typer.Argument(..., help="Run folder (with plan.json and mapping.yaml)."),
              approve: bool = typer.Option(False, "--approve"),
              public: bool = typer.Option(False, "--public"), synthetic: bool = typer.Option(False, "--synthetic"),
              confirm_deidentified: bool = typer.Option(False, "--confirm-deidentified"),
              confirm_use_permitted: bool = typer.Option(False, "--confirm-use-permitted"),
              data_use_statement: Optional[str] = typer.Option(None, "--data-use-statement"),
              irb_or_dua: Optional[str] = typer.Option(None, "--irb-or-dua")) -> None:
    """Approve (with --approve) and apply the run folder's mapping.yaml: BYOD folder, OMOP export, BYOD validation and
    the data-use check. Nothing is ingested; continue with `measure-it launch run <input> --out <run>`."""
    from .orchestrate import run_launch
    pj = run / "plan.json"
    if not pj.is_file():
        _out(f"no plan.json in {run}: create the run with `measure-it launch map <input> --out {run}`")
        raise typer.Exit(2)
    p = json.loads(pj.read_text())
    res = run_launch(p["input"], p["condition"], p["measurement"], out=run, llm="none", approve=approve,
                     public=public or p["options"].get("public", False),
                     synthetic=synthetic or p["options"].get("synthetic", False),
                     data_use=_data_use(confirm_deidentified, confirm_use_permitted, data_use_statement, irb_or_dua),
                     outreach=False, stop_after="data_use_gate", echo=None)
    _summary(res)
    st = res["plan"]["status"]
    if st.startswith("partial") or st == "completed":
        _out(f"\nnext: `measure-it launch run {p['input']} --condition {p['condition']} --measurement "
             f"{p['measurement']} --out {run}`")
    raise typer.Exit(0 if st.startswith("partial") or st == "completed" else (3 if st == "stopped" else 1))


@app.command("report")
def report_cmd(run: Path = typer.Argument(..., help="Run folder with launch.json.")) -> None:
    """Re-render LAUNCH_REPORT.md from the run folder's launch.json (deterministic template)."""
    from .report import rerender
    if not (run / "launch.json").is_file():
        _out(f"no launch.json in {run}")
        raise typer.Exit(2)
    _out(rerender(run)["report_md"])
