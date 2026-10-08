"""Pure helpers: NUCC taxonomy code -> specialty / facility groups.

The groups live in configs/specialty_groups.yaml (hand-chosen NUCC codes with
their NUCC display names). Nothing here touches the network or the NPPES file,
so every function is unit-testable.

Semantics (see the YAML header):
* an NPI matches a group if ANY of its taxonomy codes is in the group;
* it matches the group "as primary" if its primary taxonomy code is in the group;
* groups overlap on purpose, so counts must never be summed across groups.

An NPI taxonomy is self-reported. It does not mean the provider evaluates or
treats Long COVID, ME/CFS, POTS or any other invisible illness.
"""
from __future__ import annotations

import re
from collections.abc import Iterable, Sequence
from pathlib import Path

import pandas as pd
import yaml

from ..config import CONFIGS

SPECIALTY_GROUPS_PATH = CONFIGS / "specialty_groups.yaml"

# The 13 groups named in SPEC section F, in SPEC order.
SPEC_SPECIALTY_GROUPS = (
    "cardiology", "electrophysiology", "neurology", "rheumatology", "immunology_allergy",
    "infectious_disease", "pulmonology", "gastroenterology", "pain_medicine", "pmr",
    "primary_care", "gynecology", "vascular",
)
GROUP_TYPES = {"specialty", "facility", "research"}
ROLES = {"physician", "nurse_practitioner", "physician_assistant", "other_individual", "facility"}

# NUCC codes are 10 characters: 9 alphanumerics then a literal X.
_CODE_RE = re.compile(r"^[0-9A-Z]{9}X$")


def normalize_taxonomy_code(code) -> str | None:
    """Upper-case and strip a taxonomy code; return None if blank or malformed."""
    if code is None:
        return None
    try:
        if pd.isna(code):
            return None
    except (TypeError, ValueError):
        pass
    s = str(code).strip().upper()
    return s if _CODE_RE.match(s) else None


def load_specialty_groups(path: Path | str | None = None) -> dict:
    """Load and structurally validate the specialty-group config."""
    with open(path or SPECIALTY_GROUPS_PATH) as fh:
        cfg = yaml.safe_load(fh)
    problems = validate_group_config(cfg)
    if problems:
        raise ValueError("invalid specialty_groups config: " + "; ".join(problems))
    return cfg


def validate_group_config(cfg: dict) -> list[str]:
    """Structural checks that do not need the NUCC file. Returns a list of problems."""
    problems: list[str] = []
    groups = cfg.get("groups") or {}
    if not groups:
        return ["no groups defined"]
    for gid in SPEC_SPECIALTY_GROUPS:
        if gid not in groups:
            problems.append(f"SPEC group {gid!r} missing")
    for gid, g in groups.items():
        if g.get("group_type") not in GROUP_TYPES:
            problems.append(f"{gid}: bad group_type {g.get('group_type')!r}")
        codes = g.get("codes") or []
        if not codes:
            problems.append(f"{gid}: no codes")
        seen = set()
        for c in codes:
            code = c.get("code")
            if normalize_taxonomy_code(code) != code:
                problems.append(f"{gid}: malformed code {code!r}")
            if code in seen:
                problems.append(f"{gid}: duplicate code {code}")
            seen.add(code)
            if not c.get("display_name"):
                problems.append(f"{gid}: {code} has no display_name")
            if c.get("role") not in ROLES:
                problems.append(f"{gid}: {code} has bad role {c.get('role')!r}")
        for x in g.get("excluded_considered") or []:
            if x.get("code") in seen:
                problems.append(f"{gid}: {x.get('code')} is both included and excluded")
    return problems


def validate_against_nucc(cfg: dict, nucc: pd.DataFrame) -> list[str]:
    """Check every included/excluded code exists in the NUCC file and names match verbatim.

    nucc must have columns 'Code' and 'Display Name' (the raw NUCC CSV layout).
    """
    names = dict(zip(nucc["Code"], nucc["Display Name"]))
    problems = []
    for gid, g in cfg["groups"].items():
        for c in list(g.get("codes") or []) + list(g.get("excluded_considered") or []):
            code = c["code"]
            if code not in names:
                problems.append(f"{gid}: {code} not in NUCC file")
            elif names[code] != c["display_name"]:
                problems.append(f"{gid}: {code} display name {c['display_name']!r} != NUCC {names[code]!r}")
    return problems


def group_code_frame(cfg: dict) -> pd.DataFrame:
    """Long table: one row per (specialty_group, code)."""
    rows = []
    for gid, g in cfg["groups"].items():
        for c in g["codes"]:
            rows.append({
                "specialty_group": gid,
                "group_label": g.get("label", gid),
                "group_type": g["group_type"],
                "spec_listed": bool(g.get("spec_listed", False)),
                "taxonomy_code": c["code"],
                "taxonomy_display_name": c["display_name"],
                "role": c["role"],
            })
    return pd.DataFrame(rows)


def code_to_groups(cfg: dict) -> dict[str, tuple[str, ...]]:
    """Map each taxonomy code to the (sorted) groups that list it."""
    out: dict[str, set[str]] = {}
    for gid, g in cfg["groups"].items():
        for c in g["codes"]:
            out.setdefault(c["code"], set()).add(gid)
    return {k: tuple(sorted(v)) for k, v in out.items()}


def groups_for_codes(codes: Iterable, mapping: dict[str, tuple[str, ...]]) -> list[str]:
    """Sorted union of the groups for a provider's taxonomy codes (malformed/blank codes ignored)."""
    groups: set[str] = set()
    for code in codes:
        c = normalize_taxonomy_code(code)
        if c:
            groups.update(mapping.get(c, ()))
    return sorted(groups)


def select_primary_taxonomy(codes: Sequence, switches: Sequence) -> tuple[str | None, str]:
    """Pick the primary taxonomy from NPPES slots 1..15.

    Returns (code, method). method is 'primary_switch_Y' when exactly one slot is
    flagged Y, 'single_code' when no slot is flagged but only one code exists,
    'first_listed_no_primary_flag' / 'first_Y_of_multiple' otherwise, and
    'no_taxonomy' when there are no codes.
    """
    pairs = [(normalize_taxonomy_code(c), (str(s).strip().upper() if s is not None and not _isna(s) else ""))
             for c, s in zip(codes, switches)]
    pairs = [(c, s) for c, s in pairs if c]
    if not pairs:
        return None, "no_taxonomy"
    ys = [c for c, s in pairs if s == "Y"]
    if len(ys) == 1:
        return ys[0], "primary_switch_Y"
    if len(ys) > 1:
        return ys[0], "first_Y_of_multiple"
    if len({c for c, _ in pairs}) == 1:
        return pairs[0][0], "single_code"
    return pairs[0][0], "first_listed_no_primary_flag"


def _isna(x) -> bool:
    try:
        return bool(pd.isna(x))
    except (TypeError, ValueError):
        return False


def role_from_nucc(grouping: str, classification: str, section: str) -> str:
    """Derive the provider role used in the YAML from NUCC grouping/classification/section."""
    if section == "Non-Individual":
        return "facility"
    if grouping == "Allopathic & Osteopathic Physicians":
        return "physician"
    if classification == "Nurse Practitioner":
        return "nurse_practitioner"
    if classification == "Physician Assistant":
        return "physician_assistant"
    return "other_individual"


def join_groups(groups: Iterable[str]) -> str:
    """Pipe-joined, sorted, de-duplicated group string ('' when none)."""
    return "|".join(sorted(set(g for g in groups if g)))
