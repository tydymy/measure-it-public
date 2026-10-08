"""Pipeline entry point: NHANES 2011-2014 and 2003-2006 -> person_concepts__nhanes / __nhanes0306.

    uv run python -m measure_it.harmonize.run
"""
from __future__ import annotations

import json


def run(log=print) -> dict:
    from .nhanes import run as nhanes_run
    return nhanes_run(echo=log)


if __name__ == "__main__":
    print(json.dumps(run(), indent=1))
