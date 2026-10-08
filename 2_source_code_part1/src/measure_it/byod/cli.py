"""`measure-it byod ...` commands (registered on the main CLI in measure_it.cli)."""
from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Optional

import typer

app = typer.Typer(no_args_is_help=True, add_completion=False,
                  help="Bring your own data: validate, ingest, evaluate, deploy and remove a user dataset "
                       "(docs/BRING_YOUR_OWN_DATA.md).")


def _out(msg: str) -> None:
    typer.echo(msg)
    sys.stdout.flush()


@app.command("validate")
def validate_cmd(
    folder: Path = typer.Argument(..., help="Dataset folder with manifest.yaml and participants.csv."),
    min_cases: Optional[int] = typer.Option(None, "--min-cases", help="Minimum cases (default 10; never below 5)."),
    min_controls: Optional[int] = typer.Option(None, "--min-controls", help="Minimum controls (default 10)."),
    as_json: bool = typer.Option(False, "--json", help="Print the report as JSON."),
) -> None:
    """Schema, privacy and geography checks; refuses (exit 1) and lists every reason; else says what ingest creates."""
    from .checks import check_dir
    rep = check_dir(folder, min_cases, min_controls)
    _out(json.dumps(rep.to_dict(), indent=1, default=str) if as_json else rep.text())
    raise typer.Exit(0 if rep.ok else 1)


@app.command("ingest")
def ingest_cmd(
    folder: Path = typer.Argument(...),
    min_cases: Optional[int] = typer.Option(None, "--min-cases"),
    min_controls: Optional[int] = typer.Option(None, "--min-controls"),
) -> None:
    """Validate, then write namespaced person-layer partitions, the registry entry and DATA_AUDIT.md, and rebuild the
    canonical unions and this dataset's Digital Phenotype Vector."""
    from .ingest import Refused, ingest
    try:
        res = ingest(folder, min_cases=min_cases, min_controls=min_controls, echo=_out)
    except Refused as e:
        _out(e.report.text())
        raise typer.Exit(1)
    _out(f"ingested {res['dataset_id']}: next `measure-it byod evaluate {folder}`")


@app.command("evaluate")
def evaluate_cmd(
    folder: Path = typer.Argument(...),
    amend: bool = typer.Option(False, "--amend", help="Accept a changed plan as a recorded amendment."),
    jobs: Optional[int] = typer.Option(None, "--jobs", help="Parallel permutation workers (default: cores - 1, <= 8)."),
) -> None:
    """Lock the manifest's primary analysis (sha256, before any outcome), then run it: CV AUROC / AUPRC with a
    participant bootstrap, sensitivity at the pre-specified specificity, label-permutation null, covariate baseline."""
    from .evaluate import evaluate
    try:
        res = evaluate(folder, amend=amend, jobs=jobs, echo=_out)
    except RuntimeError as e:
        _out(f"byod evaluate refused: {e}")
        raise typer.Exit(1)
    _out(f"evaluated {res['dataset_id']}: next `measure-it byod deploy {folder}"
         + (" --demo`" if res["record"]["demo"] else "`"))


@app.command("subgroup")
def subgroup_cmd(
    folder: Path = typer.Argument(..., help="Dataset folder whose manifest has a subgroup_analysis section."),
    amend: bool = typer.Option(False, "--amend", help="Accept a changed subgroup plan as a recorded amendment."),
    jobs: Optional[int] = typer.Option(None, "--jobs", help="Parallel permutation workers (default: cores - 1, <= 8)."),
) -> None:
    """Lock the manifest's subgroup plan (sha256, before any outcome), then, within the label's cases, model
    device-positive vs device-negative on harmonized non-device features (L1 logistic, repeated CV, bootstrap CI,
    permutation null) and export the computable phenotype + strata fractions (aggregate only)."""
    from .subgroup import subgroup
    try:
        res = subgroup(folder, amend=amend, jobs=jobs, echo=_out)
    except RuntimeError as e:
        _out(f"byod subgroup refused: {e}")
        raise typer.Exit(1)
    ph = res["phenotype"]
    _out(f"phenotype {ph['phenotype_id']}: {len(ph['features'])} feature(s) "
         f"({ph['ehr_evaluable_summary']['n_ehr_evaluable']} EHR-evaluable), AUROC {ph['performance']['auroc']:.3f} "
         f"[{ph['performance']['decision']}]; results/byod/{res['dataset_id'][len('byod_'):]}/SUBGROUP.md")


@app.command("similar")
def similar_cmd(target: str = typer.Argument(..., help="Dataset folder or dataset id.")) -> None:
    """Which people resemble the device-positive patients (NHANES estimates + OMOP/clinic query packs; aggregate)."""
    from ..similar.cli import byod_similar
    raise typer.Exit(byod_similar(target, echo=_out))


@app.command("deploy")
def deploy_cmd(
    folder: str = typer.Argument(..., help="Dataset folder or dataset id."),
    demo: bool = typer.Option(False, "--demo", help="Include demo-tagged (synthetic) records in the rankings."),
    export_static: bool = typer.Option(False, "--export-static", help="Also write a static snapshot to "
                                                                       "results/byod/<id>/static_snapshot/."),
    build_db: bool = typer.Option(False, "--build-db", help="Also rebuild data/processed/measure_it_public.duckdb."),
) -> None:
    """Metric link + evidence_weighted re-scoring with the user record(s); report how the ranking moved."""
    from .deploy import deploy
    try:
        deploy(folder, demo=demo, export_static=export_static, build_db=build_db, echo=_out)
    except RuntimeError as e:
        _out(f"byod deploy refused: {e}")
        raise typer.Exit(1)


@app.command("remove")
def remove_cmd(
    dataset_id: str = typer.Argument(..., help="Dataset id (with or without the byod_ prefix)."),
    no_rescore: bool = typer.Option(False, "--no-rescore", help="Do not rebuild the metric link and rankings."),
    purge_ledger: bool = typer.Option(False, "--purge-ledger", help="Also delete the dataset's ledger lines."),
    build_db: bool = typer.Option(False, "--build-db", help="Also rebuild the DuckDB."),
) -> None:
    """Remove a user dataset and all its records; restore the rankings without it."""
    from .remove import remove
    try:
        res = remove(dataset_id, rescore=not no_rescore, purge_ledger=purge_ledger, build_db=build_db, echo=_out)
    except RuntimeError as e:
        _out(f"byod remove refused: {e}")
        raise typer.Exit(1)
    _out(f"removed {res['dataset_id']}; re-scored {len(res.get('rescored', []))} combination(s)")


@app.command("list")
def list_cmd(as_json: bool = typer.Option(False, "--json")) -> None:
    """User datasets on this machine (status, plan hash, record)."""
    from .tools import list_user_datasets
    env = list_user_datasets()
    if as_json:
        _out(json.dumps(env, indent=1, default=str))
        return
    rows = (env.get("data") or {}).get("datasets") or []
    if not rows:
        _out("no user datasets (measure-it byod ingest <dir>)")
        return
    for r in rows:
        _out(f"{r['dataset_id']}: {r['title']} [{r['status']}] owner {r['owner']}; {r['primary_label']} -> "
             f"{r['primary_condition_id']}; class {r['measurement_class']}; comparator {r['comparator_type']}; "
             f"demo {r['demo']}; plan {str(r.get('plan_sha256') or '-')[:12]}")


@app.command("show")
def show_cmd(dataset_id: str = typer.Argument(..., help="Dataset id (with or without the byod_ prefix).")) -> None:
    """The dataset's locked plan, models, performance record and last deploy report (get_user_dataset_metric JSON;
    aggregate only)."""
    from .tools import get_user_dataset_metric
    env = get_user_dataset_metric(dataset_id)
    _out(json.dumps(env, indent=1, default=str))
    raise typer.Exit(0 if env.get("status") == "ok" else 1)
