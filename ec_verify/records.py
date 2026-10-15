# ec_verify/records.py
"""The frozen cell table ec_verify/records.json (docs/records.md §1.3).

cells["<torsion_key>/<rank>"] = {torsion, torsion_key, rank, record, known_count, class_count,
exact_count}. `record` is the smallest certified known curve that reaches the cell at the freeze
(certified at that rank or higher: list fewer points); its height is the cell's baseline, the
reference the Yukon total uses, and it never changes. A key absent from the table is a
breakthrough cell. The verifier reads version and record only.
"""

from __future__ import annotations

import json
import pathlib
import re

from .cells import cell_key, torsion_key
from .score import SCORING_VERSION

_TABLE = None
_KEY = re.compile(r"^(1|[0-9]+|[0-9]+x[0-9]+)/[0-9]+$")


def load(path: pathlib.Path | None = None) -> dict:
    global _TABLE
    if path is None and _TABLE is not None:
        return _TABLE
    p = path or pathlib.Path(__file__).with_name("records.json")
    t = json.loads(p.read_text())
    if path is None:
        _TABLE = t
    return t


def cell(table: dict, key: str) -> dict | None:
    return table["cells"].get(key)


def baseline(table: dict, key: str) -> float | None:
    """The cell's scoring reference: its record height, or None for a cell the table does not have."""
    rec = record(table, key)
    return None if rec is None else float(rec["naive_height"])


def record(table: dict, key: str) -> dict | None:
    c = cell(table, key)
    return None if c is None else c.get("record")


def validate(table: dict) -> None:
    for k in ("version", "scoring_version", "ledger_frozen_at", "built_at", "cells"):
        if k not in table:
            raise ValueError(f"records table lacks {k}")
    if table["scoring_version"] != SCORING_VERSION:
        raise ValueError(f"scoring_version {table['scoring_version']} != {SCORING_VERSION}")
    if int(table["version"]) < 0:
        raise ValueError("version must be >= 0")
    for key, c in table["cells"].items():
        if not _KEY.match(key):
            raise ValueError(f"bad cell key {key!r}")
        if key != cell_key(torsion_key(c["torsion"]), c["rank"]) or c["torsion_key"] != torsion_key(c["torsion"]):
            raise ValueError(f"cell key {key!r} disagrees with its torsion/rank")
        rec = c.get("record")
        if not rec or not all(f in rec for f in ("id", "naive_height", "ainvs", "certified_rank")) or len(rec["ainvs"]) != 5:
            raise ValueError(f"cell {key}: record missing or malformed")
        if not float(rec["naive_height"]) > 0:
            raise ValueError(f"cell {key}: record naive_height must be positive (contribution divides by it)")
        if int(rec["certified_rank"]) < int(c["rank"]):
            raise ValueError(f"cell {key}: record certified_rank {rec['certified_rank']} does not reach rank {c['rank']}")
        if int(c.get("known_count", 0)) < 1 or int(c.get("class_count", 0)) < 1:
            raise ValueError(f"cell {key}: known_count/class_count must be >= 1")
