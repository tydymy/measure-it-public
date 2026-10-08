"""Launch agent part C: the orchestrated run (docs/AGENT_CONTRACT.md section 4).

    run_launch(input_path, condition, measurement, out=None, llm=None, approve=False, public=False, top=10,
               per_county=20, keep=False, ...) -> dict

Steps (each recorded in `<out>/plan.json` as done / skipped + why / failed + error / stopped + message, with timings):

    profile -> map -> approve -> apply -> omop_export -> byod_validate -> data_use_gate -> byod_ingest ->
    byod_evaluate -> byod_subgroup -> byod_similar -> query_pack_check -> byod_deploy -> geo_subgroup_burden ->
    launch_ranking -> outreach -> teardown -> report

What the data hold decides the steps: no label -> no evaluate; no `subgroup_analysis.device_positive` -> no subgroup,
no similar, no subgroup burden; a demo (synthetic) record -> no deploy unless `demo_rankings=True` for this run;
private data whose manifest does not confirm `data_use` -> the run STOPS before ingest and says which manifest fields
to confirm. An unapproved mapping stops the run before apply (contract rule 3).

Side-effect hygiene. `byod ingest`, `byod evaluate`, `byod subgroup` and `byod deploy` write shared processed tables.
The run records the dataset id, refuses to touch a dataset id that is already ingested (it never tears down data it
did not create) and, unless `keep=True`, ends with `byod remove <id>` (which restores the rankings when the run
deployed). Engine outputs under results/byod/<id>/ are copied into `<out>/engine_results/` before the teardown. The
device-subgroup burden is computed in memory (geography.subgroup.build) and written only under the run folder; the
shared `geo_subgroup_burden` table is never written by a launch run. The run folder carries a `.gitignore` of `*`.

Parts A (profile, mapper, llm) and B (apply, omop, omop_run) are called through the thin adapters below by the module
names of the contract; an input that already is a BYOD folder (manifest.yaml + participants.csv) skips map/apply.
"""
from __future__ import annotations

import getpass
import importlib
import inspect
import json
import re
import shutil
import sys
import time
import traceback
from pathlib import Path
from typing import Any, Callable

import yaml

from ..config import PROJECT_ROOT, utc_now_iso

RUN_ROOT = PROJECT_ROOT / "launch_runs"
STEPS = ["profile", "map", "approve", "apply", "omop_export", "byod_validate", "data_use_gate", "byod_ingest",
         "byod_evaluate", "byod_subgroup", "byod_similar", "query_pack_check", "byod_deploy", "geo_subgroup_burden",
         "launch_ranking", "outreach", "teardown", "report"]
# steps that still run after a stop (the report says why the run stopped)
AFTER_STOP = ("report",)
LLM_CHOICES = ("none", "claude")
FAST_ANALYSIS = {"n_repeats": 3, "n_permutations": 20, "n_bootstrap": 200}
PLACEHOLDER_RE = re.compile(r"\b(todo|tbd|to confirm|placeholder|fill in|xxx)\b", re.I)

# ------------------------------------------------------------------------------------------------------------------
# adapters to parts A and B (the functions of measure_it.agent.{read,profile,mapper,apply,omop,omop_run})
# ------------------------------------------------------------------------------------------------------------------
MAPPER_FUNCS = ["propose_mapping", "build_mapping", "make_mapping", "propose", "map_profile", "run"]


class AdapterMissing(RuntimeError):
    """A part A/B function is not available (module not built yet, or no known function name)."""


def _module(name: str):
    try:
        return importlib.import_module(f"measure_it.agent.{name}")
    except ImportError as e:
        raise AdapterMissing(f"measure_it.agent.{name} is not available ({e})") from e


def part_profile(input_path: Path, out: Path) -> tuple[dict, Any]:
    """Part A: read.read_input -> profile.profile_dataset -> profile.write_profile(<out>/profile.json)."""
    R, P = _module("read"), _module("profile")
    ds = R.read_input(input_path)
    prof = P.profile_dataset(ds)
    P.write_profile(prof, out)
    return prof, ds


def part_map(profile: dict, dataset, condition: str, measurement: str, llm: str, out: Path, echo=None) -> dict:
    """Part A: mapper.propose_mapping (rules, dictionary, fuzzy, optional Claude) -> mapping dict."""
    M = _module("mapper")
    f = next((getattr(M, n) for n in MAPPER_FUNCS if callable(getattr(M, n, None))), None)
    if f is None:
        raise AdapterMissing(f"measure_it.agent.mapper has none of {MAPPER_FUNCS}")
    res = call_with(f, profile=profile, prof=profile, dataset=dataset, ds=dataset, condition=condition,
                    condition_hint=condition, measurement=measurement, measurement_hint=measurement,
                    llm=None if llm == "none" else llm, provider=None if llm == "none" else llm,
                    use_llm=llm == "claude", out=out, out_dir=out, run_dir=out,
                    echo=echo or (lambda *a, **k: None))
    mapping = _as_dict(res, out / "mapping.yaml")
    w = getattr(M, "write_mapping", None)
    if mapping and callable(w):
        w(mapping, out)                     # part A's writer (with its reviewer header)
    return mapping


def part_apply(dataset, mapping: dict, out: Path, approve: bool, public: bool) -> dict:
    """Part B: apply.apply_mapping(tables, mapping, out) -> <out>/byod/, <out>/omop/, <out>/apply_report.json."""
    A = _module("apply")
    tables = {t.name: t.df for t in dataset.tables}
    return A.apply_mapping(tables, mapping, out, approve=approve, public=public)


def part_omop_from_byod(byod_dir: Path, omop_dir: Path) -> dict:
    """Part B: OMOP export of an existing BYOD folder (input that already is a BYOD folder)."""
    return _module("omop").export_omop_from_byod(byod_dir, omop_dir)


def part_query_pack(omop_dir: Path, phenotype_path: Path) -> dict:
    """Part B: run the stage-D All of Us query pack against the export -> <omop_dir>/query_pack_check.json."""
    return _module("omop_run").check_query_pack(omop_dir, phenotype_path)


PARTS = {"profile": part_profile, "mapper": part_map, "apply": part_apply, "omop": part_omop_from_byod,
         "omop_run": part_query_pack}


def resolve(part: str) -> Callable:
    """The part A/B adapter for `part` (tests replace entries of PARTS or monkeypatch this function)."""
    return PARTS[part]


def call_with(f: Callable, **aliases) -> Any:
    """Call f with the keyword arguments its signature accepts. `aliases` maps every name we can supply (several
    spellings of the same value are fine); a required parameter we cannot supply raises TypeError with its name."""
    sig = inspect.signature(f)
    params = sig.parameters
    if any(p.kind == p.VAR_KEYWORD for p in params.values()):
        named = {k: v for k, v in aliases.items() if k in params}
        extra = {k: v for k, v in aliases.items() if k not in params}
        return f(**named, **{k: v for k, v in extra.items() if not k.startswith("_")})
    kw = {}
    for name, p in params.items():
        if p.kind in (p.VAR_POSITIONAL, p.VAR_KEYWORD):
            continue
        if name in aliases:
            kw[name] = aliases[name]
        elif p.default is p.empty:
            raise TypeError(f"{getattr(f, '__module__', '?')}.{f.__name__}: cannot supply required parameter "
                            f"{name!r}")
    pos_only = [n for n, p in params.items() if p.kind == p.POSITIONAL_ONLY and n in kw]
    args = [kw.pop(n) for n in pos_only]
    return f(*args, **kw)


def _as_dict(obj, fallback: Path | None = None) -> dict:
    """A part's return value -> dict: a dict, a path to JSON/YAML, an object with to_dict(), or the fallback file."""
    if isinstance(obj, dict):
        return obj
    if hasattr(obj, "to_dict") and callable(obj.to_dict):
        d = obj.to_dict()
        if isinstance(d, dict):
            return d
    for p in (obj, fallback):
        if isinstance(p, (str, Path)) and Path(p).is_file():
            txt = Path(p).read_text()
            return json.loads(txt) if str(p).endswith(".json") else (yaml.safe_load(txt) or {})
    return {}


# ------------------------------------------------------------------------------------------------------------------
# plan
# ------------------------------------------------------------------------------------------------------------------
class Skip(Exception):
    """The step is not applicable to these data (recorded as skipped + reason)."""


class Stop(Exception):
    """The run must stop here (recorded as stopped + message; only AFTER_STOP steps run afterwards)."""


class Plan:
    def __init__(self, out: Path, header: dict):
        self.path = out / "plan.json"
        if self.path.is_file():                 # an earlier run of this folder: keep its plan under logs/
            (out / "logs").mkdir(exist_ok=True)
            shutil.move(str(self.path), out / "logs" / f"plan_{time.strftime('%Y%m%d-%H%M%S')}.json")
        self.data = {**header, "status": "running", "started_at": utc_now_iso(), "finished_at": None,
                     "steps": [], "message": None}
        self.stopped: str | None = None
        self.halted: str | None = None          # stop_after reached (a partial run on request, not a stop)
        self.failed: set[str] = set()
        self.write()

    def write(self) -> None:
        self.path.write_text(json.dumps(self.data, indent=1, default=str))

    def status(self, name: str) -> str | None:
        return next((s["status"] for s in self.data["steps"] if s["step"] == name), None)

    def done(self, name: str) -> bool:
        return self.status(name) == "done"

    def record(self, name: str, status: str, t0: float, **kw) -> None:
        self.data["steps"].append({"step": name, "status": status, "started_at": kw.pop("started_at", None),
                                   "seconds": round(time.time() - t0, 2), **kw})
        self.write()

    def run(self, name: str, fn: Callable[[], dict | None], requires: tuple = ()) -> dict | None:
        t0, started = time.time(), utc_now_iso()
        if self.stopped and name not in AFTER_STOP:
            self.record(name, "skipped", t0, started_at=started, reason=f"run stopped at {self.stopped}")
            return None
        if self.halted and name not in AFTER_STOP + ("teardown",):
            self.record(name, "skipped", t0, started_at=started, reason=f"not requested (partial run up to "
                                                                        f"{self.halted})")
            return None
        bad = [r for r in requires if not self.done(r)]
        if bad:
            why = "; ".join(f"{r} {self.status(r) or 'not run'}" for r in bad)
            self.record(name, "skipped", t0, started_at=started, reason=f"needs {', '.join(bad)} ({why})")
            return None
        try:
            out = fn() or {}
        except Skip as e:
            self.record(name, "skipped", t0, started_at=started, reason=str(e))
            return None
        except Stop as e:
            self.stopped = name
            self.data["message"] = str(e)
            self.record(name, "stopped", t0, started_at=started, message=str(e))
            return None
        except Exception as e:  # noqa: BLE001 - recorded, the run continues with what does not depend on it
            self.failed.add(name)
            tb = traceback.format_exc().strip().splitlines()[-6:]
            self.record(name, "failed", t0, started_at=started, error=f"{type(e).__name__}: {e}", traceback=tb)
            return None
        self.record(name, "done", t0, started_at=started, outputs=_jsonable(out))
        return out

    def provisional_status(self) -> str:
        st = [s["status"] for s in self.data["steps"]]
        if self.stopped:
            return "stopped"
        if "failed" in st:
            return "completed_with_failures"
        return f"partial (up to {self.halted})" if self.halted else "completed"

    def finish(self) -> None:
        self.data["status"] = self.provisional_status()
        self.data["finished_at"] = utc_now_iso()
        self.write()


def _jsonable(o):
    return json.loads(json.dumps(o, default=str))


# ------------------------------------------------------------------------------------------------------------------
# helpers
# ------------------------------------------------------------------------------------------------------------------
def local_only(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    gi = path / ".gitignore"
    if not gi.exists():
        gi.write_text("# measure-it launch: local-only run folder (may hold person-level derived data, named "
                      "clinicians); never commit\n*\n")
    return path


def default_out(condition: str, measurement: str) -> Path:
    stamp = time.strftime("%Y%m%d-%H%M%S")
    return RUN_ROOT / f"{stamp}_{condition}_{measurement}"


def is_byod_folder(p: Path) -> bool:
    return p.is_dir() and (p / "manifest.yaml").is_file() and (p / "participants.csv").is_file()


class Log:
    """Engine chatter goes to <out>/logs/engine.log (and to echo when verbose)."""

    def __init__(self, path: Path, echo: Callable[[str], None] | None):
        path.parent.mkdir(parents=True, exist_ok=True)
        self.path, self.echo = path, echo

    def __call__(self, *msg) -> None:
        s = " ".join(str(m) for m in msg)
        with open(self.path, "a") as fh:
            fh.write(s + "\n")
        if self.echo:
            self.echo(s)


def data_use_problems(m: dict) -> list[str]:
    """Manifest fields the owner must confirm before any ingest (empty list = confirmed)."""
    du = m.get("data_use") or {}
    out = []
    if du.get("deidentified") is not True:
        out.append("data_use.deidentified: true   (the data are de-identified under the rules that apply to them)")
    if du.get("use_permitted") is not True:
        out.append("data_use.use_permitted: true  (an IRB protocol / data-use agreement / consent allows this local "
                   "secondary analysis)")
    st = str(du.get("statement") or "").strip()
    if not st or PLACEHOLDER_RE.search(st):
        out.append("data_use.statement: <one sentence naming the permission>")
    return out


def read_manifest(byod_dir: Path) -> dict:
    return yaml.safe_load((byod_dir / "manifest.yaml").read_text()) or {}


def write_manifest(byod_dir: Path, m: dict, note: str) -> None:
    (byod_dir / "manifest.yaml").write_text(f"# {note}\n" + yaml.safe_dump(m, sort_keys=False, width=110,
                                                                           allow_unicode=True))


def patch_manifest(byod_dir: Path, *, synthetic: bool, public: bool, input_name: str, data_use: dict | None,
                   carried: dict | None, fast: bool, measurement: str) -> list[str]:
    """The run's own, recorded edits to the manifest part B wrote: synthetic/demo flags, data_use for synthetic or
    public inputs or from the user's explicit confirmation (or a confirmation carried over from an earlier run of the
    same folder), fast analysis settings, the measurement bundle. Returns what was changed."""
    m = read_manifest(byod_dir)
    changed = []
    du = dict(m.get("data_use") or {})
    if synthetic:
        if not (m.get("synthetic") and m.get("demo")):
            changed.append("synthetic: true, demo: true (--synthetic)")
        m["synthetic"], m["demo"] = True, True
        du.update({"deidentified": True, "use_permitted": True})
        if not str(du.get("statement") or "").strip() or PLACEHOLDER_RE.search(str(du.get("statement"))):
            du["statement"] = "SYNTHETIC: generated data, no real person; no IRB or data-use agreement applies."
        du.setdefault("irb_or_dua_reference", "not applicable (SYNTHETIC)")
    elif m.get("synthetic"):
        if not m.get("demo"):
            changed.append("demo: true (a synthetic dataset is always demo)")
        m["demo"] = True
    if public and not synthetic:
        if data_use_problems({"data_use": du}):
            changed.append("data_use for a public dataset (--public)")
        du.update({"deidentified": True, "use_permitted": True})
        if not str(du.get("statement") or "").strip() or PLACEHOLDER_RE.search(str(du.get("statement"))):
            du["statement"] = (f"Public, already de-identified dataset ({input_name}) used under its published terms "
                               "(--public; attested by the person running the launch).")
    if carried and not synthetic and not public and data_use_problems({"data_use": du}) \
            and not data_use_problems({"data_use": carried}):
        du.update(carried)
        changed.append("data_use carried over from this run folder's earlier, confirmed manifest")
    if data_use:
        upd = {k: v for k, v in data_use.items() if v is not None}
        if upd:
            du.update(upd)
            changed.append(f"data_use fields confirmed on the command line: {sorted(upd)}")
    m["data_use"] = du
    if fast:
        for sec in ("primary_analysis", "subgroup_analysis"):
            if isinstance(m.get(sec), dict):
                m[sec].update(FAST_ANALYSIS)
        changed.append(f"fast analysis settings {FAST_ANALYSIS} (--fast; smoke runs only)")
    meas = m.get("measurement") or {}
    if isinstance(meas, dict) and meas.get("class") and not meas.get("bundle") and measurement != meas.get("class"):
        try:
            from ..scoring.opportunity import measurement_spec
            spec = measurement_spec(measurement)
            if spec.get("kind") == "bundle" and meas["class"] in spec.get("members", []):
                meas["bundle"] = spec["measurement_id"]
                m["measurement"] = meas
                changed.append(f"measurement.bundle: {spec['measurement_id']} (contains class {meas['class']})")
        except Exception:  # noqa: BLE001 - leave the manifest as part B wrote it
            pass
    if changed:
        write_manifest(byod_dir, m, "manifest written by measure-it launch apply; edited by the launch run: "
                       + "; ".join(changed))
    return changed


def save_mapping(path: Path, mapping: dict) -> None:
    """Write mapping.yaml, keeping the leading comment header of the existing file (part A's reviewer notes)."""
    header = ""
    if path.is_file():
        lines = path.read_text().splitlines(keepends=True)
        header = "".join(ln for ln in lines[:next((i for i, ln in enumerate(lines) if not ln.startswith("#")),
                                                   len(lines))])
    path.write_text(header + yaml.safe_dump(mapping, sort_keys=False, width=120, allow_unicode=True))


def prepare_mapping(m: dict, *, synthetic: bool, data_use: dict | None, fast: bool) -> list[str]:
    """The run's recorded edits to the mapping before apply (part B builds the manifest from it): synthetic flag,
    data_use confirmed on the command line (or the SYNTHETIC statement), fast analysis settings."""
    edits = []
    who = f"{getpass.getuser()} at {utc_now_iso()}"
    if synthetic and not m.get("synthetic"):
        m["synthetic"] = True
        edits.append("synthetic: true (--synthetic)")
    du = dict(m.get("data_use") or {})
    if synthetic and data_use_problems({"data_use": du}):
        du.update({"deidentified": True, "use_permitted": True,
                   "statement": "SYNTHETIC: generated data, no real person; no IRB or data-use agreement applies.",
                   "irb_or_dua_reference": "not applicable (SYNTHETIC)"})
        m["data_use_confirmed_via"] = f"--synthetic ({who})"
        edits.append("data_use: SYNTHETIC statement (--synthetic)")
    upd = {k: v for k, v in (data_use or {}).items() if v is not None}
    if upd:
        du.update(upd)
        m["data_use_confirmed_via"] = f"measure-it launch --confirm-* / --data-use-statement ({who})"
        edits.append(f"data_use fields confirmed on the command line: {sorted(upd)}")
    if du:
        m["data_use"] = du
    if fast:
        m["primary_analysis"] = {**(m.get("primary_analysis") or {}), **FAST_ANALYSIS}
        if isinstance(m.get("subgroup_analysis"), dict):
            m["subgroup_analysis"].update(FAST_ANALYSIS)
        edits.append(f"fast analysis settings {FAST_ANALYSIS} (--fast; smoke runs only)")
    return edits


def _dataset_present(ds: str) -> bool:
    from ..byod import common as K
    from ..byod.person import dataset_row
    try:
        known = dataset_row(ds) is not None
    except Exception:  # noqa: BLE001
        known = False
    return known or bool(K.dataset_partitions(ds)) or K.raw_path(ds).exists()


# ------------------------------------------------------------------------------------------------------------------
# the run
# ------------------------------------------------------------------------------------------------------------------

DEVICE_POSITIVE_RE = re.compile(r"^\s*([A-Za-z0-9_.\-]+)\s*(?:(>=|<=|>|<)\s*(-?\d+(?:\.\d+)?))?\s*$")


def device_positive_plan(mapping: dict, spec: str, condition: str) -> dict:
    """`--device-positive COLUMN[>=X]` -> a pre-specified subgroup_analysis (locked by `byod subgroup` before any
    outcome). COLUMN is a source column mapped to role `device`; without a comparison it must be a 1/0 flag."""
    from .apply import safe_feature, safe_name
    m = DEVICE_POSITIVE_RE.match(spec or "")
    if not m:
        raise ValueError(f"--device-positive {spec!r}: use COLUMN (a 1/0 flag) or COLUMN>=X / >X / <=X / <X")
    col, op, thr = m.group(1), m.group(2) or ">=", float(m.group(3)) if m.group(3) else 1.0
    devices = {}
    for t in (mapping.get("tables") or {}).values():
        for c, cs in (t.get("columns") or {}).items():
            if isinstance(cs, dict) and cs.get("role") == "device":
                devices[c] = cs
    if col not in devices:
        raise ValueError(f"--device-positive: {col!r} is not a column mapped to role 'device'; device columns: "
                         f"{sorted(devices) or 'none'} (edit mapping.yaml or re-run with --remap)")
    cs = devices[col]
    block = safe_name(cs.get("block") or "device", "device")
    return {"name": f"{block}_positive", "within_label": condition,
            "description": f"{condition} cases with {col} {op} {thr:g} (device-positive), pre-specified on the command line",
            "device_positive": {"block": f"device_{block}" if not block.startswith("device_") else block,
                                "feature": safe_feature(cs.get("feature") or col), "threshold": thr, "direction": op},
            "ehr_domains": ["condition", "measurement", "drug", "demographic", "survey"], "max_features": 10}


def run_launch(input_path, condition: str, measurement: str, out=None, llm: str | None = None, approve: bool = False,
               public: bool = False, top: int = 10, per_county: int = 20, *, keep: bool = False,
               synthetic: bool = False, demo_rankings: bool = False, data_use: dict | None = None,
               fast: bool = False, outreach: bool = True, remap: bool = False, stop_after: str | None = None,
               verbose: bool = False, echo: Callable[[str], None] | None = print,
               device_positive: str | None = None) -> dict:
    """Run the launch agent end to end; returns {'out', 'plan', 'launch', 'report'}. See the module docstring."""
    input_path = Path(input_path).resolve()
    llm = (llm or "none").lower()
    if llm not in LLM_CHOICES:
        raise ValueError(f"--llm must be one of {LLM_CHOICES}")
    if keep and demo_rankings:
        raise ValueError("--demo-rankings puts a demo record into the rankings for this run only; it cannot be "
                         "combined with --keep (demo records never stay in the default rankings)")
    if not input_path.exists():
        raise FileNotFoundError(input_path)
    if stop_after is not None and stop_after not in STEPS:
        raise ValueError(f"stop_after must be one of {STEPS}")
    out = Path(out).resolve() if out else default_out(condition, measurement)
    local_only(out)
    if not (out.parent / ".gitignore").exists() and out.parent == RUN_ROOT:
        local_only(RUN_ROOT)
    say = echo or (lambda *_: None)
    log = Log(out / "logs" / "engine.log", echo if verbose else None)
    header = {"run": out.name, "out": str(out), "input": str(input_path), "condition": condition,
              "measurement": measurement,
              "options": {"llm": llm, "approve": approve, "public": public, "top": top, "per_county": per_county,
                          "keep": keep, "synthetic": synthetic, "demo_rankings": demo_rankings, "fast": fast,
                          "outreach": outreach, "remap": remap, "stop_after": stop_after},
              "dataset_id": None, "teardown": None}
    plan = Plan(out, header)
    ctx: dict[str, Any] = {"out": out, "input": input_path, "condition": condition, "measurement": measurement,
                           "llm": llm, "byod_input": is_byod_folder(input_path), "ingested": False,
                           "deployed": False, "log": log}
    prev_manifest = out / "byod" / "manifest.yaml"
    ctx["carried_data_use"] = None
    if prev_manifest.is_file() and not ctx["byod_input"]:
        try:
            ctx["carried_data_use"] = (yaml.safe_load(prev_manifest.read_text()) or {}).get("data_use")
        except yaml.YAMLError:
            pass

    def step(name, fn, requires=()):
        say(f"[launch] {name} ...")
        r = plan.run(name, fn, requires)
        if stop_after == name:
            plan.halted = name
        s = plan.data["steps"][-1]
        say(f"[launch] {name}: {s['status']}" + (f" ({s.get('reason') or s.get('error') or s.get('message')})"
                                                  if s["status"] != "done" else f" in {s['seconds']}s"))
        return r

    # ---- A: profile, map, approve
    def do_profile():
        if ctx["byod_input"]:
            raise Skip("input is already a BYOD folder (manifest.yaml + participants.csv): profiled by `byod "
                       "validate` instead")
        prof, dataset = resolve("profile")(input_path, out)
        if not prof:
            raise RuntimeError("part A profile returned nothing")
        if not (out / "profile.json").is_file():
            (out / "profile.json").write_text(json.dumps(prof, indent=1, default=str))
        ctx["profile"], ctx["dataset"] = prof, dataset
        return {"profile_json": str(out / "profile.json"), "format": prof.get("format"),
                "tables": [{"name": t.get("name"), "n_rows": t.get("n_rows"), "n_columns": len(t.get("columns", []))}
                           for t in prof.get("tables", [])]}

    def do_map():
        if ctx["byod_input"]:
            raise Skip("input is already a BYOD folder: nothing to map")
        mp = out / "mapping.yaml"
        if mp.is_file() and not remap:
            # the human-editable mapping of this run folder is never overwritten (re-run with --remap to discard it)
            old = yaml.safe_load(mp.read_text()) or {}
            if device_positive and not old.get("subgroup_analysis"):
                if old.get("approved_by"):
                    raise Stop("--device-positive would change an approved mapping.yaml; add subgroup_analysis by "
                               "hand and re-approve, or re-run with --remap")
                old["subgroup_analysis"] = device_positive_plan(old, device_positive, condition)
                save_mapping(mp, old)
            ctx["mapping"] = old
            return {"mapping_yaml": str(mp), "summary": mapping_summary(old),
                    "reused": ("existing mapping.yaml of this run folder, approved by "
                               f"{old['approved_by']} at {old.get('approved_at')}") if old.get("approved_by") else
                              "existing mapping.yaml of this run folder (not yet approved)"}
        mapping = call_with(resolve("mapper"), profile=ctx["profile"], dataset=ctx.get("dataset"),
                            condition=condition, measurement=measurement, llm=llm, out=out, echo=log)
        if not mapping:
            raise RuntimeError("part A mapper returned no mapping")
        mapping.setdefault("condition_hint", condition)
        mapping.setdefault("measurement_hint", measurement)
        mapping.setdefault("approved_by", None)
        mapping.setdefault("approved_at", None)
        if device_positive and not mapping.get("subgroup_analysis"):
            mapping["subgroup_analysis"] = device_positive_plan(mapping, device_positive, condition)
        save_mapping(mp, mapping)
        ctx["mapping"] = mapping
        return {"mapping_yaml": str(mp), "summary": mapping_summary(mapping)}

    def do_approve():
        if ctx["byod_input"]:
            raise Skip("no mapping (BYOD folder input)")
        m = ctx["mapping"]
        if m.get("approved_by"):
            return {"approved_by": m["approved_by"], "approved_at": m.get("approved_at"), "how": "already approved"}
        if not approve:
            raise Stop(f"mapping.yaml is not approved. Review {out / 'mapping.yaml'} (every column's role, code and "
                       f"the proposed_by rationale), then run `measure-it launch apply {out} --approve` or re-run "
                       "`measure-it launch run ... --approve` with the same --out.")
        m["approved_by"] = f"{getpass.getuser()} (measure-it launch --approve)"
        m["approved_at"] = utc_now_iso()
        save_mapping(out / "mapping.yaml", m)
        return {"approved_by": m["approved_by"], "approved_at": m["approved_at"], "how": "--approve"}

    # ---- B: apply, OMOP
    def do_apply():
        if ctx["byod_input"]:
            ctx["byod_dir"] = input_path
            raise Skip("input is already a BYOD folder: used as is")
        m = ctx["mapping"]
        edits = prepare_mapping(m, synthetic=synthetic, data_use=data_use, fast=fast)
        if edits:
            save_mapping(out / "mapping.yaml", m)
        res = resolve("apply")(ctx["dataset"], m, out, approve, public)
        bd = out / "byod"
        if not (bd / "manifest.yaml").is_file():
            raise RuntimeError(f"part B apply did not write {bd}/manifest.yaml")
        ctx["byod_dir"] = bd
        ctx["apply_result"] = res if isinstance(res, dict) else {}
        changes = patch_manifest(bd, synthetic=synthetic, public=public, input_name=input_path.name,
                                 data_use=data_use, carried=ctx["carried_data_use"], fast=fast,
                                 measurement=measurement)
        ar = ctx["apply_result"]
        return {"byod_dir": str(bd), "files": sorted(p.name for p in bd.iterdir() if p.is_file()),
                "mapping_edits_by_launch": edits, "manifest_edits_by_launch": changes,
                "apply_ok": ar.get("ok"), "data_use_source": ar.get("data_use_source"),
                "participant_id_transform": ar.get("participant_id_transform"),
                "n_participants": ar.get("n_participants")}

    def do_omop():
        od = out / "omop"
        if (od / "person.parquet").exists() and not ctx["byod_input"]:
            return {"omop_dir": str(od), "written_by": "apply", "tables": omop_tables(od)}
        omop_err = (ctx.get("apply_result") or {}).get("omop") or {}
        if omop_err.get("ok") is False:
            raise RuntimeError(f"OMOP export inside apply failed: {omop_err.get('error')}")
        resolve("omop")(ctx["byod_dir"], od)
        if not (od / "person.parquet").exists():
            raise RuntimeError(f"part B omop export did not write {od}/person.parquet")
        return {"omop_dir": str(od), "written_by": "omop.export_omop_from_byod", "tables": omop_tables(od)}

    # ---- BYOD validate, data-use gate
    def do_validate():
        bd = ctx.get("byod_dir")
        if bd is None:
            raise Skip("no BYOD folder")
        if ctx["byod_input"] and (synthetic or public or data_use or fast):
            raise RuntimeError("the input BYOD folder is used as is; --synthetic/--public/--confirm-*/--fast edit a "
                               "manifest only when the launch run wrote it (edit your manifest.yaml instead)")
        from ..byod.checks import check_dir
        rep = check_dir(bd)
        d = rep.to_dict()
        (out / "byod_validation.json").write_text(json.dumps(d, indent=1, default=str))
        m = read_manifest(bd)
        ctx["manifest"] = m
        ctx["ds"] = None
        if m.get("dataset_id"):
            from ..byod.common import engine_id
            ctx["ds"] = engine_id(str(m["dataset_id"]))
            plan.data["dataset_id"] = ctx["ds"]
        ctx["validation_ok"] = bool(rep.ok)
        res = {"ok": bool(rep.ok), "n_errors": len(d.get("errors", [])), "n_warnings": len(d.get("warnings", [])),
               "report": str(out / "byod_validation.json")}
        if not rep.ok:
            ctx["validation_errors"] = d.get("errors", [])
            raise RuntimeError(f"BYOD validation refused the folder ({len(d.get('errors', []))} reason(s)): "
                               + "; ".join(str(e) for e in d.get("errors", [])[:5]))
        return res

    def do_gate():
        m = ctx.get("manifest")
        if m is None:
            raise Skip("no manifest")
        probs = data_use_problems(m)
        kind = "synthetic" if m.get("synthetic") else ("public" if public else "private")
        if probs:
            where = (f"{ctx['byod_dir'] / 'manifest.yaml'}" if ctx["byod_input"] else
                     f"the data_use block of {out / 'mapping.yaml'} (part B writes the manifest from it)")
            raise Stop(f"{kind} data: the data use is not confirmed, so nothing is ingested. As the data owner, "
                       f"confirm these fields in {where}, or pass --confirm-deidentified --confirm-use-permitted "
                       "--data-use-statement '...' [--irb-or-dua REF], and re-run with the same --out: "
                       + " | ".join(probs))
        if not ctx.get("validation_ok"):
            raise Skip("BYOD validation failed (see byod_validate): nothing is ingested")
        ds = ctx.get("ds")
        if ds and _dataset_present(ds):
            raise Stop(f"{ds} is already ingested on this machine; a launch run never modifies or tears down a "
                       f"dataset it did not create. Remove it first (`measure-it byod remove {ds}`) or set another "
                       "dataset_id in mapping.yaml.")
        return {"data_kind": kind, "synthetic": bool(m.get("synthetic")), "demo": bool(m.get("demo")),
                "statement": (m.get("data_use") or {}).get("statement")}

    # ---- engine layers
    def do_ingest():
        from ..byod.ingest import Refused, ingest
        ctx["ingested"] = True                      # a partial ingest is torn down too
        plan.data["teardown"] = {"pending": not keep, "keep": keep, "dataset_id": ctx["ds"]}
        plan.write()
        try:
            res = ingest(ctx["byod_dir"], echo=log)
        except Refused as e:
            raise RuntimeError("ingest refused: " + "; ".join(e.report.errors[:5])) from e
        return {"dataset_id": res["dataset_id"], "partitions": res["partitions"]}

    def do_evaluate():
        m = ctx["manifest"]
        pa = m.get("primary_analysis") or {}
        labels = [lb.get("column") for lb in m.get("labels") or []]
        if not pa.get("label") or pa.get("label") not in labels:
            raise Skip("no label (manifest has no primary_analysis.label among its labels): nothing to evaluate")
        from ..byod.evaluate import evaluate
        res = evaluate(ctx["byod_dir"], echo=log)
        rec = res["record"]
        ctx["record"] = rec
        return {"record_id": rec.get("record_id"), "auroc": rec.get("auroc"), "auroc_ci_low": rec.get("auroc_ci_low"),
                "auroc_ci_high": rec.get("auroc_ci_high"), "decision": rec.get("decision"), "demo": rec.get("demo")}

    def do_subgroup():
        sa = ctx["manifest"].get("subgroup_analysis") or {}
        dp = sa.get("device_positive") or {}
        if not (dp.get("column") or (dp.get("block") and dp.get("feature") is not None)):
            raise Skip("no device_positive definition (manifest subgroup_analysis.device_positive): no subgroup")
        from ..byod.subgroup import subgroup
        res = subgroup(ctx["byod_dir"], echo=log)
        ph = res["phenotype"]
        ctx["phenotype"] = ph
        return {"phenotype_id": ph["phenotype_id"], "n_features": len(ph["features"]),
                "auroc": ph["performance"]["auroc"], "decision": ph["performance"]["decision"],
                "pooled_fraction": ph["strata_fractions"][0]["fraction"]}

    def do_similar():
        if "phenotype" not in ctx:
            raise Skip("no computable phenotype (no subgroup)")
        from ..similar.cli import run as similar_run
        s = similar_run(ctx["byod_dir"], log=log)
        ctx["similar"] = s
        return {"phenotypes": [e["phenotype_id"] for e in s["phenotypes"]],
                "references": {e["phenotype_id"]: {k: v.get("status") for k, v in e["references"].items()}
                               for e in s["phenotypes"]}}

    def do_query_pack_check():
        from ..byod import common as K
        od = out / "omop"
        packs = K.results_path(ctx["ds"]) / "similar"
        if not (od / "person.parquet").exists():
            raise Skip("no OMOP export")
        if not packs.is_dir() or not any(packs.iterdir()):
            raise Skip("no stage-D query pack (no computable phenotype)")
        # the pack is regenerated from the phenotype by part B (same generator as `byod similar`), then executed
        ph_paths = sorted(K.results_path(ctx["ds"]).glob("computable_phenotype*.json"))
        if not ph_paths:
            raise Skip("no computable_phenotype.json")
        res = resolve("omop_run")(od, ph_paths[0])
        qp = od / "query_pack_check.json"
        if not qp.exists() and isinstance(res, dict):
            qp.write_text(json.dumps(res, indent=1, default=str))
        if not qp.exists():
            raise RuntimeError("part B omop_run did not write omop/query_pack_check.json")
        d = json.loads(qp.read_text())
        if d.get("status") == "sqlglot_missing":
            raise Skip("sqlglot is not installed (not a project dependency): the pack was generated but not executed; "
                       "re-run with `uv run --with sqlglot measure-it launch run ...` to execute it")
        return {"query_pack_check": str(qp), "summary": query_pack_summary(d)}

    def do_deploy():
        rec = ctx.get("record")
        if not rec:
            raise Skip("no performance record (no evaluate)")
        if rec.get("demo") and not demo_rankings:
            raise Skip("demo (synthetic) record: it stays out of the rankings (MEASURE_IT_BYOD_DEMO semantics); the "
                       "launch ranking uses the engine's ranking without it. --demo-rankings includes it for this "
                       "run only (the teardown restores the rankings)")
        from ..byod.deploy import deploy
        ctx["deployed"] = True
        rep = deploy(ctx["byod_dir"], demo=bool(rec.get("demo") and demo_rankings), echo=log)
        from ..byod import common as K
        return {"demo_mode": rep.get("demo_mode"), "changed_scored_combinations": rep.get("changed_scored_combinations"),
                "deploy_report": str(K.results_path(ctx["ds"]) / "DEPLOY_REPORT.md")}

    def do_geo_subgroup():
        ph = ctx.get("phenotype")
        if not ph:
            raise Skip("no device subgroup (no computable phenotype)")
        from .launch import subgroup_burden
        sub, notes = subgroup_burden(ctx["ds"], ph["phenotype_id"])
        if sub is None or sub.empty:
            raise RuntimeError(f"no subgroup burden rows: {notes}")
        path = out / "geo_subgroup_burden.csv"
        sub.to_csv(path, index=False)
        ctx["subgroup_burden"] = sub
        nat = sub[sub["geo_level"] == "national"]
        return {"rows": int(len(sub)), "csv": str(path), "condition_id": str(sub["condition_id"].iloc[0]),
                "burden_basis": sorted(sub["burden_basis"].astype(str).unique()),
                "national": nat[["subgroup_estimate", "subgroup_ci_low", "subgroup_ci_high"]].to_dict("records"),
                "shared_table_written": False, "notes": notes}

    def do_launch():
        from .launch import launch_ranking
        res = launch_ranking(condition, measurement, top=top, subgroup=ctx.get("subgroup_burden"),
                             phenotype=ctx.get("phenotype"), log=log)
        ctx["launch"] = res
        (out / "launch_ranking.json").write_text(json.dumps(res, indent=1, default=str))
        return {"primary": res["primary"], "engine_weight_set": res["engine"].get("weight_set"),
                "engine_basis": res["engine"].get("basis"),
                "top_engine": [c["geo_name"] for c in res["engine"].get("counties", [])[:5]],
                "top_subgroup": [c["geo_name"] for c in (res.get("subgroup") or {}).get("counties", [])[:5]]}

    def do_outreach():
        if not outreach:
            raise Skip("--no-outreach")
        if "launch" not in ctx:
            raise Skip("no launch ranking")
        from .launch import outreach_workbook
        res = outreach_workbook(condition, measurement, ctx["launch"], out / "outreach", top=top,
                                per_county=per_county)
        ctx["outreach"] = res
        return {k: res.get(k) for k in ("path", "county_basis", "clinicians_listed", "with_measurement_activity",
                                        "with_experience_site_match", "capacity_basis")}

    def do_teardown():
        if not ctx["ingested"]:
            raise Skip("nothing was ingested")
        ds = ctx["ds"]
        copied = copy_engine_results(ds, out / "engine_results")
        if keep:
            plan.data["teardown"] = {"pending": False, "keep": True, "dataset_id": ds,
                                     "note": f"--keep: {ds} stays ingested; remove it with `measure-it byod remove "
                                             f"{ds}`"}
            raise Skip(f"--keep: {ds} left ingested (remove with `measure-it byod remove {ds}`); engine results "
                       f"copied ({len(copied)} files)")
        res = teardown(ds, rescore=ctx["deployed"], log=log)
        res["engine_results_copied"] = len(copied)
        plan.data["teardown"] = {"pending": False, "keep": False, "dataset_id": ds, **res}
        if not res["verified"]:
            raise RuntimeError(f"teardown incomplete: {res}")
        return res

    def do_report():
        from .report import write_report
        result = assemble(ctx, plan)
        paths = write_report(out, result)
        ctx["report_paths"] = paths
        return paths

    try:
        step("profile", do_profile)
        step("map", do_map, requires=() if ctx["byod_input"] else ("profile",))
        step("approve", do_approve)
        step("apply", do_apply, requires=() if ctx["byod_input"] else ("approve",))
        if ctx["byod_input"] and "byod_dir" not in ctx:
            ctx["byod_dir"] = input_path
        step("omop_export", do_omop, requires=() if ctx["byod_input"] else ("apply",))
        step("byod_validate", do_validate)
        step("data_use_gate", do_gate)
        step("byod_ingest", do_ingest, requires=("data_use_gate",))
        step("byod_evaluate", do_evaluate, requires=("byod_ingest",))
        step("byod_subgroup", do_subgroup, requires=("byod_ingest",))
        step("byod_similar", do_similar, requires=("byod_subgroup",))
        step("query_pack_check", do_query_pack_check, requires=("byod_similar",))
        step("byod_deploy", do_deploy, requires=("byod_evaluate",))
        step("geo_subgroup_burden", do_geo_subgroup, requires=("byod_subgroup",))
        step("launch_ranking", do_launch)
        step("outreach", do_outreach, requires=("launch_ranking",))
    finally:
        # never leave user data in the engine silently, whatever happened above
        stopped = plan.stopped
        plan.stopped = None
        step("teardown", do_teardown)
        plan.stopped = stopped
    step("report", do_report)
    plan.finish()
    return {"out": str(out), "plan": plan.data, "launch": ctx.get("launch"),
            "report": ctx.get("report_paths"), "dataset_id": ctx.get("ds")}


# ------------------------------------------------------------------------------------------------------------------
# teardown, summaries, assembly
# ------------------------------------------------------------------------------------------------------------------
def copy_engine_results(ds: str, dest: Path) -> list[str]:
    from ..byod import common as K
    src = K.results_path(ds)
    if not src.is_dir():
        return []
    local_only(dest)
    out = []
    for p in sorted(src.rglob("*")):
        if p.is_file() and p.name != ".gitignore":
            q = dest / p.relative_to(src)
            q.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(p, q)
            out.append(str(q.relative_to(dest)))
    return out


def clear_engine_caches() -> int:
    """Clear every functools cache of the loaded measure_it modules. After `byod remove` the processed tables no longer
    hold the dataset, but a long-lived process (MCP server, API, a test session) would otherwise keep serving cached
    rows that cite it (e.g. its phenotype signatures in a recommendation)."""
    import sys
    n = 0
    for name, mod in list(sys.modules.items()):
        if not name.startswith("measure_it.") or mod is None:
            continue
        for obj in list(vars(mod).values()):
            if callable(getattr(obj, "cache_clear", None)) and hasattr(obj, "cache_info"):
                obj.cache_clear()
                n += 1
    return n


def teardown(ds: str, *, rescore: bool, log=print) -> dict:
    """`byod remove <ds>` and a check that nothing of the dataset is left in the processed tables."""
    from ..byod import common as K
    from ..byod.remove import remove
    from ..store import read_table, table_exists
    try:
        res = remove(ds, rescore=rescore, echo=log)
    except RuntimeError as e:
        res = {"dataset_id": ds, "removed": [], "note": str(e)}
    clear_engine_caches()
    left = [p.name for p in K.dataset_partitions(ds)]
    rows = {}
    for t in (K.DATASETS_TABLE, K.RECORDS_TABLE, K.PHENOTYPES_TABLE, K.STRATA_TABLE):
        rows[t] = int((read_table(t)["dataset_id"].astype(str) == ds).sum()) if table_exists(t) else 0
    verified = not left and not any(rows.values()) and not K.raw_path(ds).exists() and not K.results_path(ds).exists()
    return {"removed_items": len(res.get("removed", [])), "rescored": [list(x) for x in res.get("rescored", [])],
            "partitions_left": left, "rows_left": rows, "verified": verified}


def mapping_summary(mapping: dict) -> dict:
    roles, methods, dropped, unmapped = {}, {}, [], []
    for tname, t in (mapping.get("tables") or {}).items():
        cols = (t or {}).get("columns") or {}
        pb = (t or {}).get("proposed_by") or {}
        for c, spec in cols.items():
            spec = spec or {}
            r = spec.get("role", "unmapped")
            roles[r] = roles.get(r, 0) + 1
            meth = (pb.get(c) or {}).get("method") or spec.get("method") or "unrecorded"
            methods[meth] = methods.get(meth, 0) + 1
            if r == "drop":
                dropped.append({"table": tname, "column": c, "reason": spec.get("reason", "")})
            elif r == "unmapped":
                unmapped.append({"table": tname, "column": c, "reason": spec.get("reason", "")})
    return {"n_columns": sum(roles.values()), "roles": roles, "methods": methods, "dropped": dropped,
            "unmapped": unmapped}


def omop_tables(od: Path) -> dict:
    out = {}
    for p in sorted(od.glob("*.parquet")):
        try:
            import pyarrow.parquet as pq
            out[p.stem] = int(pq.ParquetFile(p).metadata.num_rows)
        except Exception:  # noqa: BLE001
            out[p.stem] = None
    return out


def query_pack_summary(d: dict) -> dict:
    """status, per-variant per-step row counts and SQL errors of a query_pack_check.json (part B)."""
    if not isinstance(d, dict):
        return {"status": None}
    out = {"status": d.get("status"), "ok": d.get("ok"), "error": d.get("error"),
           "phenotype_id": d.get("phenotype_id"), "cdm_tables_missing": d.get("cdm_tables_missing"),
           "coverage": d.get("coverage"), "variants": {}, "n_steps_with_counts": 0, "n_errors": 0, "errors": []}
    if d.get("error"):
        out["n_errors"] += 1
        out["errors"].append(str(d["error"])[:300])
    for name, v in (d.get("variants") or {}).items():
        steps = v.get("steps") or []
        errs = [s for s in steps if s.get("error")]
        fin = v.get("final") or {}
        if fin.get("error"):
            errs.append({"step": "final", "error": fin["error"]})
        out["variants"][name] = {"ok": v.get("ok"), "steps": [[s.get("step"), s.get("rows")] for s in steps],
                                 "final_rows": len(fin.get("rows") or [])}
        out["n_steps_with_counts"] += sum(1 for s in steps if s.get("rows") is not None)
        out["n_errors"] += len(errs)
        out["errors"] += [f"{name}/{e.get('step')}: {str(e.get('error'))[:200]}" for e in errs]
    out["errors"] = out["errors"][:6]
    return out


def _llm_requests(out: Path, mapping: dict | None) -> list[dict]:
    for name in ("llm_requests.jsonl", "llm_log.jsonl", "llm_requests.json"):
        p = out / name
        if p.is_file():
            txt = p.read_text().strip()
            if not txt:
                return []
            if name.endswith(".jsonl"):
                return [json.loads(x) for x in txt.splitlines() if x.strip()]
            v = json.loads(txt)
            return v if isinstance(v, list) else [v]
    reqs = (mapping or {}).get("llm_requests")
    return list(reqs) if isinstance(reqs, list) else []


def assemble(ctx: dict, plan: Plan) -> dict:
    """Everything the report needs, from the run's own files and the engine's outputs (no model writes any of it)."""
    out: Path = ctx["out"]
    m = ctx.get("manifest") or {}
    mapping = ctx.get("mapping")
    if mapping is None and (out / "mapping.yaml").is_file():
        mapping = yaml.safe_load((out / "mapping.yaml").read_text()) or {}
    prof = ctx.get("profile")
    if prof is None and (out / "profile.json").is_file():
        prof = json.loads((out / "profile.json").read_text())
    val = json.loads((out / "byod_validation.json").read_text()) if (out / "byod_validation.json").is_file() else None
    qp = None
    if (out / "omop" / "query_pack_check.json").is_file():
        qp = json.loads((out / "omop" / "query_pack_check.json").read_text())
    omop_dir = out / "omop"
    res = {
        "generated_at": utc_now_iso(), "engine": "measure-it launch (docs/AGENT_CONTRACT.md section 4)",
        "run": {**{k: plan.data.get(k) for k in ("run", "out", "input", "condition", "measurement", "options",
                                                 "dataset_id", "message", "started_at")},
                "status": plan.provisional_status()},
        "steps": [{k: s.get(k) for k in ("step", "status", "seconds", "reason", "error", "message")}
                  for s in plan.data["steps"]],
        "input": {"profile": profile_summary(prof), "mapping": mapping_summary(mapping) if mapping else None,
                  "mapping_approved_by": (mapping or {}).get("approved_by"),
                  "mapping_approved_at": (mapping or {}).get("approved_at"),
                  "llm": ctx["llm"], "llm_requests": _llm_requests(out, mapping)},
        "harmonization": {
            "byod_validation": None if val is None else {"ok": val.get("ok"), "errors": val.get("errors", []),
                                                         "warnings": val.get("warnings", [])},
            "manifest": {"dataset_id": m.get("dataset_id"), "synthetic": bool(m.get("synthetic")),
                         "demo": bool(m.get("demo")), "measurement": m.get("measurement"),
                         "comparator": m.get("comparator"),
                         "labels": [{k: lb.get(k) for k in ("column", "condition", "label_basis")}
                                    for lb in m.get("labels") or []],
                         "subgroup_analysis": bool(m.get("subgroup_analysis"))},
            "omop": {"dir": str(omop_dir), "tables": omop_tables(omop_dir)} if omop_dir.is_dir() else None,
            "query_pack_check": None if qp is None else {"path": str(out / "omop" / "query_pack_check.json"),
                                                         **query_pack_summary(qp)},
        },
        "device_evidence": record_summary(ctx.get("record")),
        "subgroup": phenotype_summary(ctx.get("phenotype")),
        "similar": similar_summary(ctx.get("similar")),
        "launch": ctx.get("launch"),
        "outreach": ctx.get("outreach"),
        "teardown": plan.data.get("teardown"),
    }
    res["caveats"] = collect_caveats(res, ctx)
    return res


def profile_summary(prof: dict | None) -> dict | None:
    if not prof:
        return None
    return {"input": prof.get("input"), "format": prof.get("format"),
            "tables": [{"name": t.get("name"), "source_file": t.get("source_file"), "n_rows": t.get("n_rows"),
                        "n_columns": len(t.get("columns") or []), "shape": t.get("shape"),
                        "id_candidates": t.get("id_candidates")} for t in prof.get("tables") or []],
            "dictionary": prof.get("dictionary")}


def record_summary(rec: dict | None) -> dict | None:
    if not rec:
        return None
    keep = ("record_id", "dataset_id", "target_condition", "primary_class", "auroc", "auroc_ci_low", "auroc_ci_high",
            "perm_p", "decision", "op_sensitivity", "op_specificity", "n_cases", "n_controls", "label_basis",
            "comparator_kind", "quality_tier", "measurement_only", "demo", "synthetic")
    out = {k: rec.get(k) for k in keep if k in rec}
    for k in ("n_cases", "n_controls"):
        if out.get(k) is not None and out[k] == out[k]:
            out[k] = int(out[k])
    cav = rec.get("caveats")
    out["caveats"] = json.loads(cav) if isinstance(cav, str) and cav.startswith("[") else (cav or [])
    return out


def phenotype_summary(ph: dict | None) -> dict | None:
    if not ph:
        return None
    return {"phenotype_id": ph["phenotype_id"], "condition_id": ph["condition_id"],
            "subgroup_definition": ph.get("subgroup_definition"), "base_population": ph.get("base_population"),
            "features": [{k: f.get(k) for k in ("feature_key", "label", "vocabulary", "coefficient_standardized",
                                                "ehr_evaluable")} for f in ph.get("features", [])],
            "performance": ph.get("performance"), "pooled_fraction": ph["strata_fractions"][0],
            "ehr_evaluable_summary": ph.get("ehr_evaluable_summary"), "caveats": ph.get("caveats", []),
            "synthetic": ph.get("synthetic"), "demo": ph.get("demo")}


def similar_summary(s: dict | None) -> dict | None:
    if not s:
        return None
    out = {"nhanes_2021_2023": s.get("nhanes_2021_2023"), "phenotypes": []}
    for e in s.get("phenotypes", []):
        refs = {}
        for k, v in e.get("references", {}).items():
            refs[k] = {"status": v.get("status"), "reason": v.get("reason"),
                       "reference_population": v.get("reference_population"), "national": v.get("national", []),
                       "notes": v.get("notes")}
        out["phenotypes"].append({"phenotype_id": e["phenotype_id"], "query_pack": e.get("query_pack"),
                                  "references": refs})
    return out


def collect_caveats(res: dict, ctx: dict) -> list[str]:
    """Every caveat the layers produced, de-duplicated, in layer order."""
    cav: list[str] = []

    def add(x):
        if x and x not in cav:
            cav.append(str(x))
    man = res["harmonization"]["manifest"]
    if man.get("synthetic"):
        add("SYNTHETIC input: generated data, no real person. Every number in this report demonstrates the mechanics "
            "only; the dataset is demo-tagged and stays out of default rankings.")
    if ctx["llm"] == "claude" and res["input"].get("llm_requests"):
        add("A language model proposed some column mappings from metadata only; every proposal was validated "
            "deterministically and approved before use. No value, score or ranking comes from a model.")
    elif ctx["llm"] == "claude":
        add("--llm claude was requested but no model request was logged (no API key or package, or nothing left "
            "unmapped; see logs/engine.log): every mapping came from rules, the dictionary or fuzzy matching.")
    mp = res["input"].get("mapping") or {}
    if mp.get("unmapped"):
        add(f"{len(mp['unmapped'])} source column(s) stayed unmapped and never entered any layer.")
    for c in (res.get("device_evidence") or {}).get("caveats", []):
        add(c)
    comp = (res.get("device_evidence") or {}).get("comparator_kind") or ((man.get("comparator") or {}).get("type"))
    if comp == "healthy" and not any("healthy" in c.lower() for c in cav):
        add("The comparator is healthy controls: device performance is optimistic; the clinical contrast is people "
            "with other causes of the same symptoms.")
    for c in (res.get("subgroup") or {}).get("caveats", []):
        add(c)
    sim = res.get("similar") or {}
    if sim.get("phenotypes"):
        add("NHANES 'similar people' estimates describe a resemblance profile in a pre-pandemic general population; "
            "they are not Long COVID prevalence. Query packs were parse-checked here and must be run inside the "
            "All of Us / N3C enclaves under their publication rules (All of Us forbids publishing cells of 1-20).")
    for c in (res.get("launch") or {}).get("caveats", []):
        add(c)
    o = res.get("outreach") or {}
    if o.get("path"):
        add("The outreach workbook names clinicians from public NPPES/CMS files; it is written only under the run "
            "folder (local-only, never committed). Billing evidence covers Medicare fee-for-service only.")
    td = res.get("teardown") or {}
    if td.get("keep"):
        add(f"--keep: the dataset {td.get('dataset_id')} is still ingested in the engine (remove it with "
            f"`measure-it byod remove {td.get('dataset_id')}`).")
    add("Every ranked county is a candidate deployment opportunity for pilot evaluation, not a validated diagnostic "
        "pathway.")
    return cav


if __name__ == "__main__":  # pragma: no cover
    from .cli import app
    sys.exit(app())
