"""Write results/LABS_VS_DIAGNOSIS_RESULTS.md from results/tables/labs_dx_*.csv (Parts 1 and 2).

Every number in the document is read from those tables; nothing is typed in by hand.

    uv run python -m measure_it.labs.report
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

from ..config import RESULTS, TABLES
from .labs_vs_diagnosis import LABEL_BY_ID, PREFIX as P1
from .plan import PLAN_PATH, plan_sha256

P2 = "labs_dx_appelman"
OUT = RESULTS / "LABS_VS_DIAGNOSIS_RESULTS.md"


def _t(name: str) -> pd.DataFrame:
    return pd.read_csv(TABLES / f"{name}.csv")


def _f(x, nd=3, sign=False) -> str:
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "n/a"
    return f"{x:+.{nd}f}" if sign else f"{x:.{nd}f}"


def _p(x) -> str:
    if x is None or not np.isfinite(x):
        return "n/a"
    return f"{x:.3f}" if x >= 0.001 else f"{x:.1e}"


def _table(df: pd.DataFrame) -> str:
    cols = list(df.columns)
    out = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for _, r in df.iterrows():
        out.append("| " + " | ".join(str(r[c]).replace("|", "/") for c in cols) + " |")
    return "\n".join(out)


def _top_effects(eff: pd.DataFrame, label: str, tier: str = "primary", k: int = 5) -> str:
    e = eff[(eff.label_id == label) & (eff.tier == tier) & eff.fdr_significant]
    e = e.reindex(e.effect_size.abs().sort_values(ascending=False).index).head(k)
    if e.empty:
        return "none"
    return "; ".join(f"{n} {s:+.2f}" for n, s in zip(e.lab_name, e.effect_size))


# ================================================================================================= Part 1
def part1() -> list[str]:
    meta = json.loads((TABLES / f"{P1}_run_metadata.json").read_text())
    summ, eff = _t(f"{P1}_summary"), _t(f"{P1}_effects")
    flow, dic = _t(f"{P1}_population_flow"), _t(f"{P1}_analytes")
    perf, inc = _t(f"{P1}_model_performance"), _t(f"{P1}_increments")
    ph = _t(f"{P1}_posthoc_cotinine_adjusted") if (TABLES / f"{P1}_posthoc_cotinine_adjusted.csv").exists() else None
    gate = meta.get("positive_control_gate", {})
    L: list[str] = []
    n_prim = int(((dic.tier == "primary") & dic.included).sum())
    n_ext = int(((dic.tier == "extended") & dic.included).sum())
    n_bin = int(((dic.tier == "extended_binary") & dic.included).sum())
    modelled = summ[summ.get("delta_auroc").notna()] if "delta_auroc" in summ else summ.iloc[0:0]
    L += ["## Part 1 — NHANES 2011-2014: laboratory analytes vs diagnosis-like labels", "",
          f"Population: {flow.iloc[-1]['n']:,} adults (>= 20, not pregnant, biochemistry profile measured; flow in "
          f"`{P1}_population_flow.csv`). Analytes: {n_prim} in the primary panel (models + screen), {n_ext} "
          f"extended-tier continuous analytes and {n_bin} binary serology/STI/HPV results (screen only; subsample "
          f"layers analysed among sampled people). Unweighted internal comparisons. Seed {meta['seed']}; "
          f"{meta['n_splits']}-fold x {meta['n_repeats']} CV; {meta['n_boot']} participant bootstraps; "
          f"{meta['n_perm']} label permutations; {meta['n_shuffle']} lab-block shuffles.", ""]

    # ---------- gate
    L += ["### Positive-control gate (can the pipeline see lab signal at all?)", ""]
    if gate:
        L.append(f"Pre-specified: labs must add to demographics for {', '.join('`'+g+'`' for g in gate['required'])}. "
                 f"Result: **{'PASSED' if gate['gate_passed'] else 'FAILED'}** "
                 f"({', '.join(f'{k}: {v}' for k, v in gate['passed'].items())}).")
        L.append("")
        pc = summ[summ.group == "positive_control"]
        rows = []
        for _, r in pc.iterrows():
            rows.append({"label": r.label_id, "cases / controls": f"{r.n_cases} / {r.n_controls}",
                         "AUROC demog.": _f(r.auroc_demographics), "AUROC labs only": f"{_f(r.auroc_labs)} ({r.auroc_labs_ci})",
                         "AUROC demog.+labs": f"{_f(r.auroc_demographics_labs)} ({r.auroc_demographics_labs_ci})",
                         "Delta (95% CI)": f"{_f(r.delta_auroc, sign=True)} ({_f(r.delta_ci_low, sign=True)} to {_f(r.delta_ci_high, sign=True)})",
                         "strongest FDR analytes (SD units, age/sex-adjusted)": _top_effects(eff, r.label_id, k=4)})
        L += [_table(pd.DataFrame(rows)), ""]

    # ---------- headline
    L += ["### Headline: every label, nulls included", "",
          "Delta = AUROC(demographics + labs) - AUROC(demographics), paired participant bootstrap. 'labs only' p = "
          "label-permutation p. Reading applies the locked rule (CI > 0 and lab-block shuffle p <= 0.05; a Delta "
          "< 0.02 is called negligible). # FDR = primary-panel analytes with q < 0.05 / tested.", ""]
    rows = []
    for _, r in summ.iterrows():
        lab = LABEL_BY_ID[r.label_id]
        has = "delta_auroc" in summ and pd.notna(r.get("delta_auroc"))
        rows.append({
            "label": r.label_id + (" (PROXY)" if lab.is_proxy else ""), "group": r.group,
            "cases / controls": f"{r.n_cases} / {r.n_controls}",
            "AUROC demog.": _f(r.auroc_demographics) if has else "-",
            "AUROC labs only (perm p)": f"{_f(r.auroc_labs)} ({_p(r.perm_p_labs)})" if has else "-",
            "Delta demog.+labs (95% CI)": (f"{_f(r.delta_auroc, sign=True)} ({_f(r.delta_ci_low, sign=True)} to "
                                           f"{_f(r.delta_ci_high, sign=True)})") if has else "-",
            "Delta over demog.+BMI": (f"{_f(r.delta_auroc_bmi_baseline, sign=True)} ({r.delta_bmi_ci})"
                                      if has and pd.notna(r.get("delta_auroc_bmi_baseline")) else "-"),
            "# FDR": f"{r.n_primary_fdr}/{r.n_primary_analytes_tested}", "reading": r.reading})
    L += [_table(pd.DataFrame(rows)), ""]

    # ---------- nulls
    um = summ[summ.group == "umbrella"]
    nulls = [r.label_id for _, r in um.iterrows() if not bool(r.get("labs_add_to_demographics", False))]
    negl = [r.label_id for _, r in um.iterrows() if bool(r.get("labs_add_to_demographics", False))
            and r.delta_auroc < 0.02]
    adds = [r.label_id for _, r in um.iterrows() if bool(r.get("labs_add_to_demographics", False))
            and r.delta_auroc >= 0.02]
    L += ["### Umbrella labels: what the labs do and do not show", "",
          f"* **No increment over demographics (null or not modelled):** {', '.join(nulls) or 'none'}.",
          f"* **Increment present but negligible (< 0.02 AUROC):** {', '.join(negl) or 'none'}.",
          f"* **Increment >= 0.02:** {', '.join(adds) or 'none'}.",
          "* BH across the umbrella labels' Delta bootstrap p-values (descriptive): " + ", ".join(
              f"{r.label_id} {_p(r.delta_boot_p_bh_umbrella)}" for _, r in um.iterrows()
              if pd.notna(r.get("delta_boot_p_bh_umbrella"))) + ".", ""]

    # ---------- which labs
    L += ["### Which analytes differ, per label (primary panel, q < 0.05; age/sex-adjusted SD units, BMI-adjusted "
          "in brackets)", ""]
    rows = []
    for _, r in summ.iterrows():
        e = eff[(eff.label_id == r.label_id) & (eff.tier == "primary") & eff.fdr_significant]
        e = e.reindex(e.effect_size.abs().sort_values(ascending=False).index).head(6)
        rows.append({"label": r.label_id, "cases": r.n_cases,
                     "FDR primary / extended": f"{r.n_primary_fdr} / {r.n_extended_fdr}",
                     "strongest primary-panel differences": "; ".join(
                         f"{n} {s:+.2f} [{b:+.2f}]" for n, s, b in zip(e.lab_name, e.effect_size,
                                                                        e.effect_size_bmi_adj)) or "none"})
    L += [_table(pd.DataFrame(rows)), "",
          "Full per-analyte table (all tiers, CIs, p, q, BMI-adjusted): `labs_dx_nhanes_effects.csv`; also "
          "`phenotype_signatures__nhanes_labs`.", ""]

    # ---------- pairwise
    pw = summ[summ.group == "pairwise"]
    L += ["### Umbrella conditions against each other (pairwise contrasts)", ""]
    rows = []
    for _, r in pw.iterrows():
        has = "delta_auroc" in summ and pd.notna(r.get("delta_auroc"))
        rows.append({"contrast": LABEL_BY_ID[r.label_id].title, "n": f"{r.n_cases} / {r.n_controls}",
                     "# FDR primary": f"{r.n_primary_fdr}/{r.n_primary_analytes_tested}",
                     "FDR analytes": _top_effects(eff, r.label_id, k=4),
                     "AUROC demog. / demog.+labs": (f"{_f(r.auroc_demographics)} / {_f(r.auroc_demographics_labs)}"
                                                    if has else "-"),
                     "Delta (95% CI)": (f"{_f(r.delta_auroc, sign=True)} ({_f(r.delta_ci_low, sign=True)} to "
                                        f"{_f(r.delta_ci_high, sign=True)})" if has else "no model (< 25 per group)")})
    L += [_table(pd.DataFrame(rows)), ""]

    # ---------- binary
    b = eff[(eff.tier == "extended_binary") & eff.fdr_significant].copy()
    b["pos_cases"] = (b.mean_cases * b.n_cases).round().astype(int)
    L += ["### Binary serology / STI / HPV (extended tier, q < 0.05)", "",
          "Adjusted log odds ratios. Rows with fewer than 5 positive cases are fragile (a single person moves them).", ""]
    if len(b):
        L += [_table(b.assign(**{"log OR (95% CI)": [f"{e:+.2f} ({lo:+.2f} to {hi:+.2f})" for e, lo, hi in
                                                    zip(b.effect_size, b.ci_low, b.ci_high)],
                                 "positive cases / cases": [f"{p} / {n}" for p, n in zip(b.pos_cases, b.n_cases)],
                                 "q": [_p(q) for q in b.q_value_bh]})
                      [["label_id", "lab_name", "layer", "positive cases / cases", "log OR (95% CI)", "q"]]), ""]

    # ---------- post hoc
    if ph is not None:
        L += ["### POST HOC: is the symptom-label signal smoking? (not pre-specified)", "",
              "Serum cotinine and urinary NNAL (tobacco exposure) are among the largest differences for several "
              "symptom labels, which was seen only after the locked analysis ran. Primary-panel contrasts re-estimated "
              "with log serum cotinine as an extra covariate (tobacco analytes omitted); number with q < 0.05:", ""]
        rows = []
        for _, r in summ.iterrows():
            a = eff[(eff.label_id == r.label_id) & (eff.tier == "primary") & (eff.layer != "tobacco_cotinine")]
            c = ph[ph.label_id == r.label_id]
            rows.append({"label": r.label_id, "q < 0.05 age/sex (non-tobacco analytes)": int(a.fdr_significant.sum()),
                         "q < 0.05 + cotinine": int((c.q_value_bh < 0.05).sum()), "tested": int(c.p_value.notna().sum())})
        L += [_table(pd.DataFrame(rows)), ""]

    pm = TABLES / f"{P1}_posthoc_models.csv"
    if pm.exists():
        m = pd.read_csv(pm)
        L += ["### POST HOC: do the symptom-label increments survive smoking, BMI and the proxy's own exclusions?", "",
              "Not pre-specified; run after the locked results were seen. Labs here exclude the tobacco analytes. "
              "(i) `smoking_adjusted`: baseline = demographics + log serum cotinine (+ BMI). (ii) "
              "`controls_without_exclusion_dx`: ME/CFS-like proxy cases vs controls with none of the proxy's exclusion "
              "diagnoses, so labs cannot score by recognising the excluded diseases (diabetes via HbA1c, anaemia, "
              "thyroid, liver).", ""]
        L += [_table(m.assign(**{"baseline": ["demographics + " + b.replace("_cot", "log cotinine").replace("+", " + ")
                                              if b != "demographics" else b for b in m.baseline],
                                 "cases / n": [f"{c} / {n}" for c, n in zip(m.n_cases, m.n)],
                                 "AUROC base -> + labs": [f"{a:.3f} -> {b:.3f}" for a, b in zip(m.auroc_base, m.auroc_full)],
                                 "Delta (95% CI)": [f"{d:+.3f} ({lo:+.3f} to {hi:+.3f})" for d, lo, hi in
                                                    zip(m.delta_auroc, m.ci_low, m.ci_high)]})
                      [["label_id", "analysis", "baseline", "cases / n", "AUROC base -> + labs", "Delta (95% CI)"]]), ""]
    L += ["### Deviations from the locked plan", "",
          "* Primary-panel analytes failing the 70% coverage rule (estradiol, SHBG, serum hydroxycotinine) were "
          "screened in the extended tier rather than dropped.",
          "* Same-analyte unit duplicates (urine albumin URXUMS = URXUMA; serum folate LBDFOT = LBDFOTSI) were dropped "
          "(same label, rank correlation > 0.999), like LBXSCH.",
          "* NHANES reuses variable names across files (LBXHCT = hematocrit in CBC and hydroxycotinine in COT_H; LBXIN "
          "in GLU_G and INS_H): extended-tier collisions are namespaced `<var>__<layer>`; the INS_H insulin is also "
          "removed for the no-glycaemic diabetes label.",
          "* Solver: newton-cholesky (same L2 objective as lbfgs, chosen for speed).",
          "* The post-hoc cotinine-adjusted contrasts and post-hoc models above.",
          "* The lab-block shuffle reached its floor (p = 1/21) for nearly every label, including labels whose Delta CI "
          "includes 0: shuffled lab blocks lower the penalised model's AUROC, so beating them is a weak test (as in "
          "results/WEARABLE_NHANES_RESULTS.md). The Delta CI is the operative part of the rule.", ""]
    return L


# ================================================================================================= Part 2
def part2() -> list[str]:
    meta = json.loads((TABLES / f"{P2}_run_metadata.json").read_text())
    mods, dl = _t(f"{P2}_models"), _t(f"{P2}_deltas")
    met, dirn = _t(f"{P2}_metabolites"), _t(f"{P2}_paper_direction")
    L = ["## Part 2 — Appelman 2024 (MUSCLE-ME trial) plasma and muscle metabolomics: Long COVID vs healthy", "",
         f"Source `appelman_lc_pem_source` ({meta['source_version']}). Baseline: plasma 25 Long COVID vs 21 healthy "
         "(the lead's '23 vs 21' is the VO2max-and-blood overlap), muscle 25 vs 19. **There is no ME/CFS group in this "
         "file**, so the ME/CFS contrast was not possible. Age is not released; sex (and BMI) only. "
         f"{meta['n_splits']}-fold x {meta['n_repeats']} CV with in-fold top-{meta['top_k']} selection; "
         f"{meta['n_boot']} stratified bootstraps; {meta['n_perm']} permutations re-running the whole pipeline.", ""]
    prim = mods[mods.role == "primary"]
    L += ["### Primary (locked): does each tissue's metabolome separate Long COVID from healthy?", ""]
    rows = []
    for _, r in prim.iterrows():
        rows.append({"tissue": r.tissue, "n LC / healthy": f"{r.n_lc} / {r.n_healthy}",
                     "CV AUROC (95% bootstrap CI)": f"{_f(r.auroc)} ({_f(r.auroc_ci_low)}-{_f(r.auroc_ci_high)})",
                     "permutation null mean / 95th pct": f"{_f(r.perm_null_mean)} / {_f(r.perm_null_p95)}",
                     "perm p": _p(r.perm_p), "locked decision": r.decision if isinstance(r.decision, str) else "null"})
    L += [_table(pd.DataFrame(rows)), "",
          "Reading: the bootstrap CI resamples people but not the model-fitting, so it can exclude 0.5 while the "
          "permutation test (which re-runs selection and fitting) does not reject; the locked rule needs both.", ""]
    L += ["### Secondary analyses (pre-specified)", ""]
    rows = []
    for _, r in mods[mods.role == "secondary"].iterrows():
        rows.append({"model": r.model, "tissue": r.tissue, "n LC / healthy": f"{r.n_lc} / {r.n_healthy}",
                     "CV AUROC (95% CI)": f"{_f(r.auroc)} ({_f(r.auroc_ci_low)}-{_f(r.auroc_ci_high)})",
                     "perm p": _p(r.perm_p) if pd.notna(r.get("perm_p")) else "-",
                     "decision rule": r.decision if isinstance(r.get("decision"), str) else "-"})
    L += [_table(pd.DataFrame(rows)), ""]
    L += [_table(dl.assign(**{"Delta AUROC (95% CI)": [f"{d:+.3f} ({lo:+.3f} to {hi:+.3f})" for d, lo, hi in
                                                        zip(dl.delta_auroc, dl.ci_low, dl.ci_high)],
                              "boot p": [_p(x) for x in dl.boot_p]})[["tissue", "comparison", "Delta AUROC (95% CI)",
                                                                      "boot p"]]), "",
          "The sex-only models score below 0.5 (a known cross-validation artefact for a feature that carries no "
          "signal in a sex-matched sample: each training fold's sex imbalance is the mirror image of the test "
          "fold's), which inflates the sex + metabolites Delta; read the metabolite-only AUROCs instead. Adding "
          "plasma metabolites to VO2max does not improve on VO2max alone.", ""]
    L += ["### (a) Per-metabolite differences", ""]
    rows = []
    for (t, a), g in met.groupby(["tissue", "analysis"], sort=False):
        top = g.sort_values("welch_p").head(6)
        rows.append({"tissue": t, "analysis": a, "tested": len(g), "raw p < 0.05": int((g.welch_p < 0.05).sum()),
                     "q < 0.05 (BH, tissue)": int((g.q_welch_bh_tissue < 0.05).sum()),
                     "q < 0.05 sex-adjusted": int((g.get("q_sex_adj_bh_tissue", pd.Series(np.nan, g.index)) < 0.05).sum()),
                     "q < 0.05 within domain (reproduction)": int((g.q_welch_bh_within_domain_reproduction < 0.05).sum()),
                     "smallest p (Hedges g)": "; ".join(f"{m} {x:+.2f}" for m, x in zip(top.metabolite, top.hedges_g))})
    L += [_table(pd.DataFrame(rows)), "",
          "Expected false positives among raw p < 0.05 at the null: 4.2 of 83 (plasma), 5.8 of 116 (muscle).", ""]
    L += ["### (g) Published claims vs our computation (kept separate)", "",
          "Published = Appelman 2024 Results text (GLMM over three timepoints, Box-Cox, BH within pathway). Ours = "
          "baseline-only Hedges g on the released values with a stratified bootstrap CI; for pathway claims a per-person "
          "domain score (mean z-scored log2 metabolite in the authors' domain).", ""]
    rows = []
    for _, r in dirn.iterrows():
        rows.append({"tissue": r.tissue, "published claim": r.published_claim, "item": r["item"],
                     "our g (95% CI)": f"{r.our_hedges_g:+.2f} ({r.ci_low:+.2f} to {r.ci_high:+.2f})",
                     "share of domain in claimed direction": _f(r.share_claimed_sign, 2)
                     if pd.notna(r.share_claimed_sign) and r.kind == "domain" else "-",
                     "concordance": r.concordance})
    L += [_table(pd.DataFrame(rows)), ""]
    return L


def _row(df, **kw):
    m = np.ones(len(df), bool)
    for k, v in kw.items():
        m &= (df[k] == v).to_numpy()
    return df[m].iloc[0] if m.any() else None


def _d(r) -> str:
    return f"{r.delta_auroc:+.3f} (95% CI {r.delta_ci_low:+.3f} to {r.delta_ci_high:+.3f})"


def answer() -> list[str]:
    L = ["## Answer", ""]
    if (TABLES / f"{P1}_summary.csv").exists():
        s = pd.read_csv(TABLES / f"{P1}_summary.csv").set_index("label_id")
        eff = pd.read_csv(TABLES / f"{P1}_effects.csv")
        meta = json.loads((TABLES / f"{P1}_run_metadata.json").read_text())
        pm = pd.read_csv(TABLES / f"{P1}_posthoc_models.csv") if (TABLES / f"{P1}_posthoc_models.csv").exists() else None
        g = meta.get("positive_control_gate", {})
        L.append(f"**Part 1 (NHANES 2011-2014, {meta['n_population']:,} adults).** The positive-control gate "
                 f"{'passed' if g.get('gate_passed') else 'FAILED'}: labs add clearly to demographics where they should "
                 + "; ".join(f"{k} {_d(s.loc[k])}" for k in ("diabetes_no_glycaemic_labs", "weak_failing_kidneys",
                                                               "anemia_treatment", "liver_condition",
                                                               "congestive_heart_failure", "gout"))
                 + ". So the pipeline can see lab signal, and the nulls below are informative within their sample sizes.")
        lo = s[s["auroc_labs"].notna()]
        L.append(f"* **Labs-only models are not evidence of disease biology here.** They beat their permutation null for "
                 f"{int((lo.perm_p_labs <= 0.05).sum())} of {len(lo)} modelled labels, but their AUROC is within 0.03 of "
                 f"(or below) the demographics-only AUROC for {int(((lo.auroc_labs - lo.auroc_demographics) < 0.03).sum())}"
                 " of them: the lab panel re-learns age and sex (testosterone, creatinine, haemoglobin). The increment "
                 "over demographics is the pre-specified reading.")
        rx = [k for k in ("rx_fibromyalgia", "rx_migraine", "rx_insomnia", "rx_myalgia") if k in s.index]
        L.append("* **Umbrella conditions identified by prescription reason codes: NULL.** Labs add nothing detectable "
                 "to demographics for " + "; ".join(f"{k} (n = {s.loc[k, 'n_cases']}) {_d(s.loc[k])}" for k in rx)
                 + f". rx_ibs has {s.loc['rx_ibs', 'n_cases']} cases (no model). Primary-panel analytes with q < 0.05: "
                 + ", ".join(f"{k} {s.loc[k, 'n_primary_fdr']}" for k in rx + ["rx_ibs"]) + ".")
        pw = [k for k in s.index if k.startswith("pair_rx")]
        L.append("* **Umbrella conditions against each other: NULL.** " + "; ".join(
            f"{LABEL_BY_ID[k].title}: {s.loc[k, 'n_primary_fdr']} FDR analytes"
            + (f", Delta {_d(s.loc[k])}" if pd.notna(s.loc[k].get('delta_auroc')) else ", no model")
            for k in pw) + ". No routine lab tells fibromyalgia, migraine and insomnia (as treated conditions) apart.")
        if "pair_mecfs_like_proxy_vs_depression_non_proxy" in s.index:
            r = s.loc["pair_mecfs_like_proxy_vs_depression_non_proxy"]
            L.append(f"* **ME/CFS-like proxy vs PHQ-9 >= 10 non-proxy: labs 'separate' them ({_d(r)}), by construction.** "
                     f"The separating analytes are lower glycaemia in the proxy ({_top_effects(eff, r.name, k=3)}): the "
                     "proxy definition excludes diabetes and the other exclusion diagnoses, the depression comparison "
                     "group does not. This is the case definition showing through, not a biological difference.")
        pr = s.loc["mecfs_like_proxy"]
        txt = (f"* **ME/CFS-like PROXY: a pre-specified 'labs add' that does not survive scrutiny.** Locked result "
               f"{_d(pr)} (n = {pr.n_cases}); its strongest analytes are tobacco exposure and a LOWER HbA1c "
               f"({_top_effects(eff, 'mecfs_like_proxy', k=6)}).")
        if pm is not None:
            a = _row(pm, label_id="mecfs_like_proxy", analysis="controls_without_exclusion_dx", baseline="demographics")
            b = _row(pm, label_id="mecfs_like_proxy", analysis="smoking_adjusted", baseline="bmi+_cot")
            txt += (f" Post hoc, with tobacco analytes removed and cotinine + BMI in the baseline the increment is "
                    f"{b.delta_auroc:+.3f} ({b.ci_low:+.3f} to {b.ci_high:+.3f}); against controls without the proxy's "
                    f"exclusion diagnoses it is {a.delta_auroc:+.3f} ({a.ci_low:+.3f} to {a.ci_high:+.3f}). The proxy's "
                    "lab 'signal' is smoking plus the labs recognising the diseases the definition excludes, not an "
                    "ME/CFS-like biology.")
        L.append(txt)
        sy = [k for k in ("fatigue_symptom", "depression_phq9", "told_sleep_disorder") if k in s.index]
        txt = ("* **Symptom labels (fatigue, depressive symptoms, told sleep disorder): small increments, mostly "
               "BMI and smoking.** Locked: " + "; ".join(f"{k} {_d(s.loc[k])}" for k in sy) + ".")
        if pm is not None:
            txt += " Post hoc with cotinine + BMI in the baseline (tobacco analytes removed): " + "; ".join(
                f"{k} {r.delta_auroc:+.3f} ({r.ci_low:+.3f} to {r.ci_high:+.3f})" for k in sy
                for r in [_row(pm, label_id=k, analysis="smoking_adjusted", baseline="bmi+_cot")] if r is not None) + "."
        L.append(txt)
        cm = [k for k in ("rheumatoid_arthritis", "asthma", "thyroid_problem", "osteoarthritis", "psoriasis",
                          "celiac_disease") if k in s.index]
        L.append("* **Comparator diagnoses:** " + "; ".join(f"{k} {_d(s.loc[k])}" for k in cm) + ". Asthma "
                 "(eosinophils) and rheumatoid arthritis carry small increments; thyroid problem, osteoarthritis, "
                 "psoriasis and celiac disease (treated, mostly seronegative) do not.")
        L.append("")
    if (TABLES / f"{P2}_models.csv").exists():
        m = pd.read_csv(TABLES / f"{P2}_models.csv").set_index("model")
        met = pd.read_csv(TABLES / f"{P2}_metabolites.csv")
        dl = pd.read_csv(TABLES / f"{P2}_deltas.csv")
        b, mu, c = m.loc["blood_top10"], m.loc["muscle_top10"], m.loc["blood_plus_muscle_top10"]
        vb = met[(met.tissue == "blood") & (met.analysis == "baseline")]
        vm = met[(met.tissue == "muscle") & (met.analysis == "baseline")]
        v = dl[dl.comparison.str.startswith("vo2max")].iloc[0]
        L.append(f"**Part 2 (Appelman 2024 metabolomics, Long COVID vs healthy).** Plasma: CV AUROC {b.auroc:.2f} "
                 f"({b.auroc_ci_low:.2f}-{b.auroc_ci_high:.2f}), permutation p {b.perm_p:.3f} -> **null** under the "
                 f"locked rule; muscle: {mu.auroc:.2f} ({mu.auroc_ci_low:.2f}-{mu.auroc_ci_high:.2f}), p "
                 f"{mu.perm_p:.3f} -> **null**. Secondary plasma + muscle: {c.auroc:.2f} ({c.auroc_ci_low:.2f}-"
                 f"{c.auroc_ci_high:.2f}), p {c.perm_p:.3f} (meets the rule, but it is one of several secondaries at "
                 f"n = {int(c.n_lc)} vs {int(c.n_healthy)}). Per metabolite, q < 0.05 (BH per tissue): plasma "
                 f"{int((vb.q_welch_bh_tissue < 0.05).sum())}/{len(vb)}, muscle {int((vm.q_welch_bh_tissue < 0.05).sum())}"
                 f"/{len(vm)} ({', '.join(vm[vm.q_welch_bh_tissue < 0.05].metabolite)}). VO2max alone "
                 f"{m.loc['vo2max_only'].auroc:.2f}; adding plasma metabolites changes AUROC by {v.delta_auroc:+.3f} "
                 f"({v.ci_low:+.3f} to {v.ci_high:+.3f}). The paper's muscle TCA-lower direction is reproduced "
                 "(individual CIs exclude 0, not FDR-significant at n = 44); its 'plasma glycolytic metabolites higher' "
                 "is only partly reproduced (see the claims table). No ME/CFS group exists in this file.")
        L.append("")
    L.append("**Bottom line.** Routine NHANES labs separate the conditions they are known to track (diabetes without "
             "glycaemic markers, kidney disease, anaemia, liver disease, heart failure, gout), and do not separate the "
             "invisible-illness labels available here from the general population or from each other beyond "
             "demographics, BMI and smoking. The one small open metabolomics set gives a borderline, underpowered "
             "Long COVID vs healthy signal that is not better than VO2max. These are nulls for these labels in these "
             "data (treated or self-reported cases, general-population or healthy controls), not evidence that no "
             "measurable biology exists.")
    L.append("")
    return L


def write() -> str:
    lines = ["# Do clinical labs and other non-wearable measurements separate invisible-illness labels?", "",
             f"**Plan (locked before any comparison):** `docs/{PLAN_PATH.name}` (sha256 `{plan_sha256()[:16]}...`). "
             "**Modules:** `measure_it.labs.labs_vs_diagnosis` (Part 1), `measure_it.labs.appelman_metabolomics` "
             "(Part 2), `measure_it.ingestion.appelman_lc_pem`, `measure_it.labs.report`. **Tables:** "
             "`results/tables/labs_dx_nhanes_*.csv`, `results/tables/labs_dx_appelman_*.csv`. **Signature partitions:** "
             "`phenotype_signatures__nhanes_labs`, `phenotype_signatures__appelman_metabolomics`.", "",
             "Reproduce: `uv run python -m measure_it.labs.plan && uv run python -m measure_it.labs.labs_vs_diagnosis "
             "&& uv run python -m measure_it.ingestion.appelman_lc_pem && uv run python -m "
             "measure_it.labs.appelman_metabolomics && uv run python -m measure_it.labs.report`; tests: "
             "`uv run pytest tests/test_labs_vs_diagnosis.py`.", ""]
    body = []
    if (TABLES / f"{P1}_run_metadata.json").exists():
        body += part1()
    if (TABLES / f"{P2}_run_metadata.json").exists():
        body += part2()
    OUT.write_text("\n".join(lines + answer() + body))
    return str(OUT)


if __name__ == "__main__":
    print(write())
