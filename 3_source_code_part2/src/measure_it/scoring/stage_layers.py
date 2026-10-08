"""Sensitivity S7-S9: the data-layer choices adopted 2026-10-07 against the rules they replaced.

Since 2026-10-07 the primary ranking uses (configs/scoring.yaml) the BRFSS 2023 MRP small-area estimate as Long COVID
burden (`long_covid_burden: brfss_sae`) and measurement-specific clinician activity as clinic_capacity where the
measurement has a dedicated billing code (`clinic_capacity_basis: activity_where_dedicated`). The variants put the
previous rules back:

S7 legacy_hps_inherited  Long COVID county burden = the Household Pulse state value inherited by every county.
S8 legacy_density        clinic_capacity = implementer-group provider and facility density (identical to the primary
                         for measurements without a dedicated code, e.g. nailfold capillaroscopy).
S9 legacy_both           S7 + S8 = the ranking before 2026-10-07.

The SAE's within-state differences come from county composition and covariates only (results/BRFSS_SAE_RESULTS.md);
its rows carry estimate_kind 'modeled_small_area', never observed county prevalence. A device-subgroup burden (geography.subgroup) is not scored here:
it scales one condition's burden by stratum fractions, and until a non-synthetic phenotype exists it would only
re-rank on demographic composition.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from ..config import RESULTS, TABLES
from ..store import read_table
from . import opportunity as O
from .report import md_table
from .sensitivity import ranking_agreement, top_set

LEGACY_MEASURE = "lc_current_pct_all_adults"        # HPS latest state value, inherited by the state's counties
COMBOS = [("long_covid_or_me_cfs", "wearable_autonomic_activity_monitoring"),
          ("long_covid_or_me_cfs", "nailfold_capillaroscopy"),
          ("long_covid", "nailfold_capillaroscopy"),
          ("long_covid", "wearable_autonomic_activity_monitoring")]
VARIANTS = {"S7_legacy_hps_inherited": (True, False), "S8_legacy_density": (False, True), "S9_legacy_both": (True, True)}
TOP_N = 10


def legacy_burden(geo: pd.DataFrame, level: str = "county") -> np.ndarray:
    b = read_table("geo_condition_burden")
    b = b[(b["condition_id"] == "long_covid") & (b["geo_level"] == level) & (b["source_measure_id"] == LEGACY_MEASURE)
          & (b["source_table"] == "geo_condition_burden__cdc_long_covid")]
    if level == "county":
        b = b[b["derivation"] == "inherited_from_state"]
    else:
        b = b[b["derivation"] == "as_published"].sort_values("period_end").groupby("geo_id").tail(1)
    return b.set_index("geo_id")["value"].reindex(geo["geo_id"]).to_numpy(float)


def variant_values(cb: O.Combo, legacy_burden_: bool) -> list:
    vals = []
    for m in cb.members:
        if legacy_burden_ and m.condition_id == "long_covid":
            vals.append(legacy_burden(cb.geo, cb.level))
        else:
            vals.append(m.value)
    return vals


def run_combo(cid: str, mid: str, level: str = "county") -> tuple[list[dict], pd.DataFrame]:
    cb = O.make_combo(cid, mid, level)
    base = O.evaluate(cb)["composite"][0]
    rank_base = O.rank_rows(base, cb.eligible)[0]
    has_lc = any(m.condition_id == "long_covid" for m in cb.members)
    rows, regions = [], []
    names = cb.geo["geo_name"].to_numpy()
    for name, (use_legacy_burden, use_density) in VARIANTS.items():
        if use_legacy_burden and not has_lc:
            continue
        cap = cb.cr["clinic_capacity_density_pct"] if use_density else None
        comp = O.evaluate(cb, member_values=variant_values(cb, use_legacy_burden), capacity_pct=cap)["composite"][0]
        elig = cb.eligible & ~np.isnan(comp)
        rank_v = O.rank_rows(comp, elig)[0]
        ag = ranking_agreement(base, comp, elig)
        t_base, t_new = set(top_set(base, elig, TOP_N)), set(top_set(comp, elig, TOP_N))
        rows.append({"analysis": name, "condition_id": cid, "measurement_id": mid, "geo_level": level,
                     "primary_capacity_basis": cb.cr.get("clinic_capacity_basis", ""), **ag,
                     "entering_top10": "; ".join(names[i] for i in sorted(t_new - t_base, key=lambda i: rank_v[i])),
                     "leaving_top10": "; ".join(names[i] for i in sorted(t_base - t_new, key=lambda i: rank_base[i]))})
        top = np.flatnonzero(((rank_base <= TOP_N) | (rank_v <= TOP_N)) & elig)
        regions.append(pd.DataFrame({"analysis": name, "condition_id": cid, "measurement_id": mid,
                                     "geo_id": cb.geo["geo_id"].to_numpy()[top], "geo_name": names[top],
                                     "rank_primary": rank_base[top], "rank_variant": rank_v[top],
                                     "composite_primary": base[top], "composite_variant": comp[top]})
                       .sort_values("rank_variant"))
    return rows, pd.concat(regions, ignore_index=True) if regions else pd.DataFrame()


def run() -> pd.DataFrame:
    rows, regs = [], []
    for cid, mid in COMBOS:
        r, g = run_combo(cid, mid)
        rows += r
        regs.append(g)
    out = pd.DataFrame(rows)
    TABLES.mkdir(parents=True, exist_ok=True)
    out.to_csv(TABLES / "stage_layers_sensitivity.csv", index=False)
    reg = pd.concat(regs, ignore_index=True)
    reg.to_csv(TABLES / "stage_layers_top10_regions.csv", index=False)
    write_report(out, reg)
    return out


def write_report(out: pd.DataFrame, reg: pd.DataFrame) -> None:
    keep = [c for c in out.columns if c not in ("entering_top10", "leaving_top10", "geo_level")]
    lines = ["# Data-layer choices of 2026-10-07 vs the rules they replaced (sensitivity S7-S9)", "",
             "_Generated by `measure_it.scoring.stage_layers`. Primary (since 2026-10-07): Long COVID burden from the "
             "BRFSS 2023 small-area model (results/BRFSS_SAE_RESULTS.md) and clinic capacity from measurement-specific "
             "Medicare billing where the measurement has a dedicated code (results/tables/stage_f_capacity_*.csv). "
             "Each variant puts a previous rule back._", "",
             "* **S7** Long COVID county burden = the inherited Household Pulse state value (no within-state "
             "information).",
             "* **S8** clinic capacity = implementer-group provider and facility density. Nailfold capillaroscopy has "
             "no billing code, so its primary already uses density and S8 changes nothing for it.",
             "* **S9** both = the ranking before 2026-10-07.", "",
             "Why the switch: the BRFSS model is newer (2023 vs HPS ending 2024-09 but state-only), reproduces CDC's "
             "published state estimates (52/52 case counts identical) and beats the national baseline in leave-one-"
             "state-out prediction (MAE 0.75 vs 0.92 pp), but its within-state ordering is not confirmed by metro-area "
             "estimates (r = 0.16). Activity capacity responds to the location shuffle (density did not, Spearman "
             "0.94-0.97) and counts clinicians who bill the measurement; it covers Medicare fee-for-service only.", "",
             "## Agreement with the primary (equal-weight) ranking", "",
             md_table(out, keep), "", "## Top-10 changes", ""]
    for r in out.itertuples():
        lines += [f"* **{r.analysis}** {r.condition_id} x {r.measurement_id}: entering {r.entering_top10 or '-'}; "
                  f"leaving {r.leaving_top10 or '-'}"]
    (RESULTS / "STAGE_LAYERS_SCORING.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    run()
