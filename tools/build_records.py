#!/usr/bin/env python3
"""Step 4 of the record-store pipeline (docs/records.md §1.3): ec_verify/records.json, the frozen cell table.

    python tools/build_records.py [--ledger ec_verify/known_curves.json] [--out ec_verify/records.json]
                                  [--previous ec_verify/records.json]

record(T, r) = the smallest certified ledger curve of torsion T with certified rank >= r. A curve
certified at rank R reaches every rank below R (list fewer of its points), so a cell exists for
every rank from 1 up to the group's top certified rank, and its record may be a higher-rank
curve (`record.certified_rank` > `rank`; `exact_count` is how many members sit at exactly r).
The record's height is the cell's baseline. The table is built once at the freeze; a rebuild is
a maintainer decision that changes every board's score, never a bot's.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

from ec_verify import cells, known_curves, records as records_mod  # noqa: E402
from ec_verify.score import SCORING_VERSION  # noqa: E402
from fetch_known import _SOURCE_ORDER  # noqa: E402

LEDGER = ROOT / "ec_verify" / "known_curves.json"
OUT = ROOT / "ec_verify" / "records.json"



def _record_key(c: dict):
    """Smallest naive height; a tie (the same curve mirrored from several sources) goes to the
    source listed first above, then to the numerically smallest id (icarm-9 before icarm-10)."""
    return (float(c["naive_height"]), _SOURCE_ORDER.get(c["source"], 9), known_curves._id_key(c["id"]))


def build_table(ledger: dict, previous: dict | None, now: str) -> dict:
    version = int(previous["version"]) + 1 if previous else 1
    prev_cells = (previous or {}).get("cells", {})
    by_torsion: dict = {}
    for c in ledger["curves"]:
        v = c["verification"]
        if v["status"] == "certified" and int(v["certified_rank"]) >= cells.MIN_SCORED_RANK:
            by_torsion.setdefault(c["torsion_key"], []).append(c)     # rank 0 proves nothing: never a cell

    def crank(c):
        return int(c["verification"]["certified_rank"])

    out: dict = {}
    for tkey, curves in by_torsion.items():
        for r in range(cells.MIN_SCORED_RANK, max(map(crank, curves)) + 1):
            members = [c for c in curves if crank(c) >= r]            # every curve that reaches rank r
            rec = min(members, key=_record_key)
            out[cells.cell_key(tkey, r)] = {
                "torsion": list(rec["torsion"]), "torsion_key": tkey, "rank": r,
                "record": {"id": rec["id"], "naive_height": float(rec["naive_height"]), "ainvs": list(rec["ainvs"]),
                           "certified_rank": crank(rec)},
                "known_count": len(members), "class_count": len({tuple(c["fingerprint"]) for c in members}),
                "exact_count": sum(1 for c in members if crank(c) == r)}
    missing = sorted(set(prev_cells) - set(out))
    if missing:
        raise SystemExit(f"!! cells present in the previous table have no certified curve now: {missing}")
    out = dict(sorted(out.items()))
    table = {"version": version, "scoring_version": SCORING_VERSION, "ledger_frozen_at": ledger["frozen_at"],
             "built_at": now, "cells": out}
    records_mod.validate(table)
    return table


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ledger", type=pathlib.Path, default=LEDGER)
    ap.add_argument("--out", type=pathlib.Path, default=OUT)
    ap.add_argument("--previous", type=pathlib.Path, default=OUT,
                    help="table whose version this build follows (default: the shipped table, whatever --out is)")
    args = ap.parse_args()
    ledger = json.loads(args.ledger.read_text())
    previous = json.loads(args.previous.read_text()) if args.previous.exists() else None
    now = dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")
    table = build_table(ledger, previous, now)
    args.out.write_text(json.dumps(table, indent=1, sort_keys=True) + "\n")
    print(f"wrote {args.out}: version {table['version']}, {len(table['cells'])} cells")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
