"""Stage D entry point: "which people resemble the device-positive patients?" (aggregate only).

    uv run python -m measure_it.similar <dataset_id | folder> [--reference nhanes --reference nhanes0306]
                                                              [--no-reference] [--no-query-packs]

For every computable phenotype of the dataset (results/byod/<id>/computable_phenotype*.json, written by
`measure-it byod subgroup`):
  1. NHANES reference cohorts: survey-weighted prevalence of the device-positive-like profile (national and by age
     band x sex), feature profile, Gower distance quantiles -> tables similar_reference_prevalence /
     similar_reference_profile / similar_reference_distance, and results/byod/<id>/similar/reference_summary__*.json.
  2. Query packs (All of Us BigQuery, N3C Spark SQL, clinic export) -> results/byod/<id>/similar/.

`byod_similar(target, echo)` is the function the `measure-it byod similar` command calls.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path


from ..config import PROCESSED
from .omop import write_query_pack
from .phenotype_io import annotate_platform, load_phenotype, phenotype_paths, resolve_target
from .reference import NHANES_2021_2023_FINDING, REFERENCES, run_reference, write_tables

DEFAULT_REFERENCES = ("nhanes", "nhanes0306")


def _local_only(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    gi = path / ".gitignore"
    if not gi.exists():
        gi.write_text("# measure-it similar: local-only material derived from a user dataset; never commit\n*\n")
    return path


def run(target: str | Path, *, references=DEFAULT_REFERENCES, reference: bool = True, query_packs: bool = True,
        log=print) -> dict:
    mid, rdir = resolve_target(target)
    paths = phenotype_paths(rdir)
    if not paths:
        raise FileNotFoundError(f"no computable_phenotype*.json under {rdir}; run `measure-it byod subgroup` first")
    out = _local_only(rdir / "similar")
    summary = {"dataset_id": mid, "phenotypes": [], "nhanes_2021_2023": NHANES_2021_2023_FINDING}
    results, ids = [], []
    for p in paths:
        ph = annotate_platform(load_phenotype(p))
        sid = re.sub(r"[^A-Za-z0-9_.-]+", "_", ph["phenotype_id"].split(":", 1)[-1])
        entry = {"phenotype_id": ph["phenotype_id"], "source": str(p), "references": {}}
        if query_packs:
            m = write_query_pack(ph, out)
            entry["query_pack"] = {k: v.get("coverage") for k, v in m["coverage"].items()}
            log(f"[similar] query pack for {ph['phenotype_id']}: " + ", ".join(
                f"{k} coverage {v['coverage']:.2f}" for k, v in m["coverage"].items()))
        if reference:
            for ds in references:
                if not (PROCESSED / f"person_concepts__{ds}.parquet").exists():
                    entry["references"][ds] = {"status": "skipped", "reason": f"person_concepts__{ds} not built "
                                                                             "(stage B+C harmonize step)"}
                    log(f"[similar] skip {ds}: person_concepts__{ds} missing")
                    continue
                res = run_reference(ph, ds, log=log)
                results.append(res)
                nat = res["prevalence"].query("age_band == 'all' and sex == 'all'")
                entry["references"][ds] = {
                    "status": "ok", "reference_population": REFERENCES[ds].label,
                    "national": nat[["variant", "prevalence", "ci_low", "ci_high", "n_unweighted",
                                     "n_positive_unweighted", "nchs_reliability", "feature_coverage", "coverage_flag",
                                     "coverage_by_domain", "weight_variable"]].to_dict("records"),
                    "notes": res["notes"]}
            ids.append(ph["phenotype_id"])
        (out / f"reference_summary__{sid}.json").write_text(json.dumps(entry, indent=1, default=str))
        summary["phenotypes"].append(entry)
    if results:
        write_tables(results, ids)
    return summary


def run_all(log=print) -> dict:
    """Pipeline entry point: stage D for every local user dataset that has a computable phenotype (no-op if none)."""
    from .phenotype_io import RESULTS_ROOT
    out = {}
    for d in sorted(p for p in RESULTS_ROOT.glob("*") if p.is_dir() and phenotype_paths(p)):
        out[d.name] = run(d.name, log=log)
    if not out:
        log("[similar] no computable phenotypes under results/byod/; nothing to do")
    return out


def byod_similar(target: str, echo=print) -> int:
    """Body of `measure-it byod similar <folder or dataset_id>`; returns the exit code."""
    try:
        s = run(target, log=echo)
    except (FileNotFoundError, KeyError, ValueError) as e:
        echo(f"byod similar refused: {e}")
        return 1
    for e in s["phenotypes"]:
        for ds, r in e["references"].items():
            if r.get("status") == "ok":
                for row in r["national"]:
                    echo(f"{e['phenotype_id']} x {ds} [{row['variant']}]: {row['prevalence']:.4f} "
                         f"({row['ci_low']:.4f}-{row['ci_high']:.4f}), coverage {row['feature_coverage']:.2f} "
                         f"({row['coverage_flag']}), {row['nchs_reliability']}")
    echo("aggregate tables: similar_reference_prevalence / _profile / _distance; files under results/byod/"
         f"{resolve_target(target)[0]}/similar/ (local-only). NHANES rows describe a resemblance profile in a "
         "pre-pandemic general population, not Long COVID prevalence.")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("target", help="dataset id (with or without byod_) or the dataset folder")
    ap.add_argument("--reference", action="append", choices=sorted(REFERENCES), help="reference cohort(s)")
    ap.add_argument("--no-reference", action="store_true", help="only write the query packs")
    ap.add_argument("--no-query-packs", action="store_true", help="only run the NHANES reference estimates")
    a = ap.parse_args(argv)
    try:
        s = run(a.target, references=tuple(a.reference or DEFAULT_REFERENCES), reference=not a.no_reference,
                query_packs=not a.no_query_packs)
    except (FileNotFoundError, KeyError, ValueError) as e:
        print(f"similar refused: {e}")
        return 1
    print(json.dumps(s, indent=1, default=str)[:4000])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
