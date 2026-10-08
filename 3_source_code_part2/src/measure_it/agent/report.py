"""Launch agent part C: `LAUNCH_REPORT.md` + `launch.json` (docs/AGENT_CONTRACT.md section 4).

A deterministic template over the structured result assembled by `orchestrate.assemble` (the run's own files and the
engine's outputs). No language model writes any text here, and every number comes from the data and the engine.
`write_report(out, result)` writes both files; `render(result)` returns the markdown; `rerender(run_dir)` rebuilds the
markdown from an existing launch.json.
"""
from __future__ import annotations

import json
from pathlib import Path

GUARDRAIL = ("Every ranked county is a candidate deployment opportunity for pilot evaluation, not a validated "
             "diagnostic pathway. Geographic joins are ecological; no person-level row is ever joined to a place.")


def _f(v, nd: int = 3) -> str:
    if v is None:
        return "-"
    if isinstance(v, bool):
        return "yes" if v else "no"
    if isinstance(v, int):
        return f"{v:,}"
    if isinstance(v, float):
        if v != v:
            return "-"
        if abs(v) >= 1000:
            return f"{v:,.0f}"
        return f"{v:.{nd}f}"
    return str(v).replace("|", "/").replace("\n", " ")


def _ci(lo, hi, nd: int = 3) -> str:
    return "-" if lo is None or hi is None else f"{_f(lo, nd)}-{_f(hi, nd)}"


def table(rows: list[dict], cols: list[tuple[str, str]]) -> list[str]:
    """cols = [(header, key or callable)]"""
    if not rows:
        return ["_none_"]
    out = ["| " + " | ".join(h for h, _ in cols) + " |", "|" + "---|" * len(cols)]
    for r in rows:
        cells = []
        for _, k in cols:
            v = k(r) if callable(k) else r.get(k)
            cells.append(v if isinstance(v, str) and callable(k) else _f(v))
        out.append("| " + " | ".join(cells) + " |")
    return out


def _counts(d: dict | None) -> str:
    return ", ".join(f"{k} {v}" for k, v in sorted((d or {}).items(), key=lambda kv: (-kv[1], kv[0]))) or "-"


def render(r: dict) -> str:
    run = r["run"]
    L: list[str] = []
    L += [f"# Launch report: {run['condition']} x {run['measurement']}", "",
          f"_Generated {r['generated_at']} by `measure-it launch` (deterministic template; no language model wrote "
          f"any text or number). Run folder `{run['out']}` (local-only). Status: **{run['status']}**._", ""]
    if run.get("message"):
        L += [f"> **The run stopped:** {run['message']}", ""]
    L += [f"> {GUARDRAIL}", ""]

    # ---- steps
    L += ["## Steps", ""]
    L += table(r["steps"], [("step", "step"), ("status", "status"), ("seconds", lambda s: _f(s.get("seconds"), 1)),
                            ("why / error", lambda s: _f(s.get("reason") or s.get("error") or s.get("message") or ""))])
    L += [""]

    # ---- what went in
    inp = r["input"]
    L += ["## 1. What went in", ""]
    p = inp.get("profile")
    if p:
        L += [f"Input `{run['input']}`, format **{p.get('format')}**"
              + (f", data dictionary {p['dictionary'].get('source')} ({p['dictionary'].get('n_fields')} fields)"
                 if isinstance(p.get("dictionary"), dict) and p["dictionary"].get("source") else "") + ".", ""]
        L += table(p.get("tables") or [], [("table", "name"), ("source file", "source_file"), ("rows", "n_rows"),
                                          ("columns", "n_columns"), ("shape", "shape"),
                                          ("id candidates", lambda t: ", ".join(t.get("id_candidates") or []) or "-")])
        L += [""]
    else:
        L += [f"Input `{run['input']}` (no profile: see the steps table).", ""]
    mp = inp.get("mapping")
    if mp:
        L += [f"Mapping: {mp['n_columns']} column(s). Roles: {_counts(mp['roles'])}. Proposed by: "
              f"{_counts(mp['methods'])}. Approved by {inp.get('mapping_approved_by') or '-'} at "
              f"{inp.get('mapping_approved_at') or '-'}.", ""]
        if mp["dropped"]:
            L += ["Dropped (never written):", ""]
            L += table(mp["dropped"], [("table", "table"), ("column", "column"), ("reason", "reason")])
            L += [""]
        if mp["unmapped"]:
            L += ["Unmapped (never entered any layer):", ""]
            L += table(mp["unmapped"], [("table", "table"), ("column", "column"), ("reason", "reason")])
            L += [""]
    if inp["llm"] == "none":
        L += ["Language model: none (`--llm none`); no model was contacted.", ""]
    else:
        reqs = inp.get("llm_requests") or []
        L += [f"Language model: `{inp['llm']}`; {len(reqs)} request(s) logged (metadata only; proposals validated "
              "deterministically before approval).", ""]

    # ---- harmonization
    h = r["harmonization"]
    L += ["## 2. Harmonization", ""]
    v = h.get("byod_validation")
    if v is not None:
        L += [f"BYOD validation (`measure-it byod validate`): **{'passed' if v['ok'] else 'refused'}**; "
              f"{len(v['errors'])} error(s), {len(v['warnings'])} warning(s)."]
        L += [f"* {e}" for e in v["errors"][:20]] + [f"* warning: {w}" for w in v["warnings"][:10]] + [""]
    man = h.get("manifest") or {}
    if man.get("dataset_id"):
        meas = man.get("measurement") or {}
        L += [f"Dataset `{man['dataset_id']}`; measurement class `{meas.get('class')}`"
              + (f" (bundle `{meas.get('bundle')}`)" if meas.get("bundle") else "")
              + f"; comparator {(man.get('comparator') or {}).get('type', '-')}; synthetic {_f(man.get('synthetic'))}, "
                f"demo {_f(man.get('demo'))}. Labels: "
              + (", ".join(f"`{lb['column']}` -> {lb['condition']} ({lb['label_basis']})" for lb in man["labels"])
                 or "none") + ".", ""]
    om = h.get("omop")
    if om:
        L += [f"All of Us-shaped OMOP CDM v5.4 export `{om['dir']}`: "
              + (", ".join(f"{k} {_f(n)}" for k, n in om["tables"].items()) or "no tables") + ".", ""]
    else:
        L += ["OMOP export: not written (see the steps table).", ""]
    qp = h.get("query_pack_check")
    if qp and qp.get("status") == "sqlglot_missing":
        L += [f"Stage-D query pack (`{qp['path']}`): generated for {qp.get('phenotype_id') or 'the phenotype'} but "
              "NOT executed: sqlglot is not installed (not a project dependency; run with `uv run --with sqlglot`).",
              ""]
    elif qp:
        L += [f"Stage-D query pack run locally against the export (`{qp['path']}`): status **{qp.get('status')}**; "
              f"{qp['n_steps_with_counts']} step(s) with row counts, {qp['n_errors']} SQL error(s); CDM tables "
              f"missing: {', '.join(qp.get('cdm_tables_missing') or []) or 'none'}."]
        for name, v in (qp.get("variants") or {}).items():
            L += [f"* {name} (ok {_f(v.get('ok'))}): " + ", ".join(f"{s_} {_f(n_)}" for s_, n_ in v.get("steps", []))
                  + f"; final rows {v.get('final_rows')}"]
        L += [f"* error: {e}" for e in qp.get("errors", [])] + [""]
    else:
        L += ["Query-pack check: not run (see the steps table).", ""]

    # ---- device evidence
    L += ["## 3. The device evidence (performance record)", ""]
    d = r.get("device_evidence")
    if d:
        L += [f"Record `{d.get('record_id')}` ({d.get('primary_class')} for {d.get('target_condition')}, tier "
              f"{_f(d.get('quality_tier'), 1)}): AUROC {_f(d.get('auroc'))} (95% CI "
              f"{_ci(d.get('auroc_ci_low'), d.get('auroc_ci_high'))}), permutation p {_f(d.get('perm_p'), 4)}, "
              f"decision **{d.get('decision')}**; sensitivity {_f(d.get('op_sensitivity'))} at specificity "
              f"{_f(d.get('op_specificity'))}; {_f(d.get('n_cases'))} cases vs {_f(d.get('n_controls'))} controls; "
              f"label basis {d.get('label_basis')}; comparator {d.get('comparator_kind')}; demo {_f(d.get('demo'))}.",
              ""]
    else:
        L += ["No performance record (no label, or evaluate not run).", ""]
    eng = (r.get("launch") or {}).get("engine") or {}
    perf = eng.get("performance") or {}
    if perf:
        L += [f"In the rankings: performance status for {run['condition']} x {run['measurement']} is "
              f"**{perf.get('status')}**" + (f" (record {perf.get('record_id')}, measurement_evidence "
                                            f"{_f(perf.get('measurement_evidence'))})"
                                            if perf.get("status") == "known" else "") + ".", ""]

    # ---- subgroup + similar
    L += ["## 4. The device subgroup and similar people", ""]
    s = r.get("subgroup")
    if s:
        perf_s = s.get("performance") or {}
        pf = s.get("pooled_fraction") or {}
        L += [f"Phenotype `{s['phenotype_id']}` (within {s['condition_id']}): {s.get('subgroup_definition')}. "
              f"Device-positive fraction among cases {_f(pf.get('fraction'))} (95% CI "
              f"{_ci(pf.get('ci_low'), pf.get('ci_high'))}, n = {_f(pf.get('n_cases'))}). Computable phenotype "
              f"(no device features): AUROC {_f(perf_s.get('auroc'))} (95% CI "
              f"{_ci(perf_s.get('auroc_ci_low'), perf_s.get('auroc_ci_high'))}), decision "
              f"**{perf_s.get('decision')}**.", ""]
        L += table(s.get("features") or [], [("feature", "feature_key"), ("label", "label"),
                                             ("vocabulary", "vocabulary"),
                                             ("std. coefficient", "coefficient_standardized"),
                                             ("EHR-evaluable", "ehr_evaluable")])
        L += [""]
    else:
        L += ["No device subgroup (no `device_positive` definition, or the subgroup step did not run).", ""]
    sim = r.get("similar")
    if sim:
        for e in sim["phenotypes"]:
            L += [f"Similar people for `{e['phenotype_id']}`; query-pack coverage: "
                  + (", ".join(f"{k} {_f(v)}" for k, v in (e.get("query_pack") or {}).items()) or "-") + "."]
            rows = []
            for ref, x in e["references"].items():
                if x.get("status") != "ok":
                    rows.append({"reference": ref, "variant": x.get("status"), "note": x.get("reason")})
                    continue
                for n in x.get("national") or []:
                    rows.append({"reference": ref, "variant": n.get("variant"), "prevalence": n.get("prevalence"),
                                 "ci": _ci(n.get("ci_low"), n.get("ci_high"), 4),
                                 "coverage": n.get("feature_coverage"), "flag": n.get("coverage_flag"),
                                 "reliability": n.get("nchs_reliability")})
            L += [""] + table(rows, [("reference", "reference"), ("variant", "variant"),
                                     ("weighted prevalence", "prevalence"), ("95% CI", "ci"),
                                     ("feature coverage", "coverage"), ("coverage flag", "flag"),
                                     ("NCHS reliability", "reliability"), ("note", "note")]) + [""]

    # ---- launch candidates
    la = r.get("launch")
    L += ["## 5. Launch candidates", ""]
    if la:
        L += [f"**{la['definition']}.**", ""]
        sub = la.get("subgroup") or {}
        if eng:
            L += [f"### (a) Engine ranking: {eng.get('basis')}"
                  + (f", weight set `{eng['weight_set']}`" if eng.get("weight_set") else " (BURDEN-ONLY)"), "",
                  f"{eng.get('why')}. " + (f"{eng.get('clinic_capacity_note')}. " if eng.get("clinic_capacity_note")
                                           else "") + f"Interval: {eng.get('monte_carlo')}.", ""]
            if eng.get("kind") == "deployment_composite":
                L += table(eng["counties"], [
                    ("rank", "rank"), ("county", "geo_name"), ("state", "state_abbr"), ("FIPS", "geo_id"),
                    ("burden", "burden_value"), ("burden pct", "burden_pct"), ("vulnerability", "vulnerability_pct"),
                    ("desert", "diagnostic_desert_pct"), ("clinic capacity", "clinic_capacity_pct"),
                    ("readiness", "research_readiness_pct"), ("device evidence", "measurement_evidence"),
                    ("composite", "composite"),
                    ("MC rank p05-p95", lambda c: _ci(c.get("mc_rank_p05"), c.get("mc_rank_p95"), 0)),
                    ("P(top 10)", "mc_p_top10")])
            else:
                keys = [k for k in (eng["counties"][0] if eng.get("counties") else {}) if k.startswith("burden_")]
                L += table(eng.get("counties") or [], [("rank", "rank"), ("county", "geo_name"),
                                                       ("state", "state_abbr"), ("FIPS", "geo_id"),
                                                       ("population", "population_total")]
                           + [(k, k) for k in keys])
            L += [""]
        if sub.get("status") == "done":
            L += [f"### (b) Device-subgroup re-score: {sub.get('basis')}"
                  + (f", weight set `{sub['weight_set']}`" if sub.get("weight_set") else " (BURDEN-ONLY)"), ""]
            if sub.get("kind") == "deployment_composite":
                eff = sub.get("effect") or {}
                L += [f"Burden member(s) replaced: {', '.join(sub['replaced_members'])} -> {sub['burden_replacement']} "
                      f"({', '.join(sub['burden_basis'])}); kept: {', '.join(sub['kept_members']) or 'none'}. Effect vs "
                      f"the engine ranking under the same weights: Spearman {_f(eff.get('spearman_vs_engine'))}, top-"
                      f"{eff.get('top_n')} overlap {eff.get('top_overlap')}; entering: "
                      f"{'; '.join(eff.get('entering') or []) or '-'}; leaving: "
                      f"{'; '.join(eff.get('leaving') or []) or '-'}. Interval: {sub.get('monte_carlo')}.", ""]
                L += table(sub["counties"], [
                    ("rank", "rank"), ("engine rank", "rank_engine_same_weights"), ("county", "geo_name"),
                    ("state", "state_abbr"), ("FIPS", "geo_id"),
                    ("subgroup adults", lambda c: _f(c.get("subgroup_estimate"), 0)),
                    ("95% CI", lambda c: _ci(c.get("subgroup_ci_low"), c.get("subgroup_ci_high"), 0)),
                    ("per 100k", lambda c: _f(c.get("subgroup_rate_per_100k"), 0)), ("burden pct", "burden_pct"),
                    ("vulnerability", "vulnerability_pct"), ("desert", "diagnostic_desert_pct"),
                    ("clinic capacity", "clinic_capacity_pct"), ("readiness", "research_readiness_pct"),
                    ("composite", "composite"),
                    ("MC rank p05-p95", lambda c: _ci(c.get("mc_rank_p05"), c.get("mc_rank_p95"), 0)),
                    ("P(top N)", "mc_p_top10")])
            else:
                L += [f"{sub.get('why')}.", ""]
                L += table(sub["counties"], [
                    ("rank", "rank"), ("county", "geo_name"), ("state", "state_abbr"), ("FIPS", "geo_id"),
                    ("subgroup cases", lambda c: _f(c.get("subgroup_estimate"), 0)),
                    ("95% CI", lambda c: _ci(c.get("subgroup_ci_low"), c.get("subgroup_ci_high"), 0)),
                    ("per 100k", lambda c: _f(c.get("subgroup_rate_per_100k"), 0))])
            L += [""]
        elif sub:
            L += [f"Device-subgroup re-score: not computed ({sub.get('reason')}).", ""]
    else:
        L += ["No launch ranking (see the steps table).", ""]

    # ---- who to contact
    L += ["## 6. Who to contact", ""]
    if la:
        prim = la["subgroup"] if la["primary"] == "subgroup" else la["engine"]
        L += table(prim.get("counties") or [], [
            ("rank", "rank"), ("county", "geo_name"), ("state", "state_abbr"),
            ("billing clinicians (Medicare FFS)", lambda c: _f((c.get("contacts") or {}).get("billing_clinicians_n"))),
            ("within 50 km", lambda c: _f((c.get("contacts") or {}).get("billing_clinicians_within_50km_n"))),
            ("code basis", lambda c: _f((c.get("contacts") or {}).get("capacity_basis"))),
            ("experience sites", lambda c: _f((c.get("contacts") or {}).get("experience_sites_n"))),
            ("candidate facilities", lambda c: _f((c.get("contacts") or {}).get("candidate_facilities_n")))])
        L += [""]
    o = r.get("outreach")
    if o and o.get("path"):
        L += [f"Outreach workbook (names, addresses, NPIs; local-only, never commit): `{o['path']}`; "
              f"{_f(o.get('clinicians_listed'))} clinician(s) listed, {_f(o.get('with_measurement_activity'))} with "
              f"measurement billing activity, {_f(o.get('with_experience_site_match'))} matched to an experience "
              f"site; county basis {o.get('county_basis')}; capacity basis {o.get('capacity_basis')}.", ""]
    else:
        L += ["Outreach workbook: not written (see the steps table).", ""]

    # ---- teardown
    td = r.get("teardown")
    L += ["## 7. Engine state after the run", ""]
    if not td:
        L += ["Nothing was ingested; the engine's processed tables were not modified.", ""]
    elif td.get("keep"):
        L += [f"`--keep`: `{td.get('dataset_id')}` is still ingested (`measure-it byod remove {td.get('dataset_id')}` "
              "removes it and restores the rankings).", ""]
    else:
        L += [f"`{td.get('dataset_id')}` was removed (`byod remove`): verified {_f(td.get('verified'))}; partitions "
              f"left {len(td.get('partitions_left') or [])}; rows left {td.get('rows_left')}; combinations re-scored "
              f"{len(td.get('rescored') or [])}. Engine results were copied to `engine_results/` in the run folder.",
              ""]

    # ---- caveats
    L += ["## 8. Caveats", ""] + [f"* {c}" for c in r.get("caveats", [])] + [""]
    return "\n".join(L)


def write_report(out: Path, result: dict) -> dict:
    out = Path(out)
    (out / "launch.json").write_text(json.dumps(result, indent=1, default=str))
    (out / "LAUNCH_REPORT.md").write_text(render(result))
    return {"launch_json": str(out / "launch.json"), "report_md": str(out / "LAUNCH_REPORT.md")}


def rerender(run_dir: Path) -> dict:
    run_dir = Path(run_dir)
    result = json.loads((run_dir / "launch.json").read_text())
    (run_dir / "LAUNCH_REPORT.md").write_text(render(result))
    return {"launch_json": str(run_dir / "launch.json"), "report_md": str(run_dir / "LAUNCH_REPORT.md")}
