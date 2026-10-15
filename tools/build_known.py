#!/usr/bin/env python3
"""Step 3 of the record-store pipeline (docs/records.md §1.2, §2): the ledger ec_verify/known_curves.json.

    python tools/build_known.py [--raw tools/known_raw.json] [--out ec_verify/known_curves.json] [--jobs 8]

Per raw entry: fingerprint (|a_p|, p <= 1000), torsion structure from gp (checked against the
source's when present), naive height and the verification summary. Fingerprint and torsion are
reused from the previous ledger when id and a-invariants are unchanged, so a rebuild only
computes the new curves. The ledger is frozen (docs/records.md §3): every entry feeds the novelty gate.
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import json
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

from ec_verify import cells, fingerprint as fp_mod, known_curves, torsion as torsion_mod  # noqa: E402
from ec_verify.curve import naive_height  # noqa: E402
from fetch_known import RAW, _SOURCE_ORDER, dump_per_line, load_raw  # noqa: E402

OUT = ROOT / "ec_verify" / "known_curves.json"
_VER_KEYS = ("status", "certified_rank", "code", "verified_at", "pari_version", "verifier_version")


def _compute(entry: dict, gp: str) -> tuple:
    fp = fp_mod.fingerprint(entry["ainvs"], gp_bin=gp)
    _, structure, _, _ = torsion_mod.strip_torsion(entry["ainvs"], [], gp_bin=gp)
    return fp, structure


def build_ledger(raw: dict, previous: dict | None, gp: str, jobs: int) -> dict:
    prev = {c["id"]: c for c in (previous or {}).get("curves", []) if c.get("torsion") is not None}
    todo, cached = [], {}
    for c in raw["curves"]:
        p = prev.get(c["id"])
        if p is not None and [int(a) for a in p["ainvs"]] == list(c["ainvs"]):
            cached[c["id"]] = (p["fingerprint"], p["torsion"])
        else:
            todo.append(c)
    with cf.ThreadPoolExecutor(max(1, jobs)) as ex:
        for c, (fp, st) in zip(todo, ex.map(lambda c: _compute(c, gp), todo)):
            cached[c["id"]] = (fp, st)
    mismatches, curves = [], []
    for c in raw["curves"]:
        fp, structure = cached[c["id"]]
        if c.get("torsion_source") is not None and list(c["torsion_source"]) != list(structure):
            mismatches.append(f"{c['id']}: source torsion {c['torsion_source']} vs gp {structure}")
        v = c.get("verification") or {"status": "unverified"}
        ver = {k: v.get(k) for k in _VER_KEYS}
        ver["status"] = v.get("status", "unverified")
        curves.append({"id": c["id"], "source": c["source"], "source_id": c["source_id"],
                       "ainvs": [str(a) for a in c["ainvs"]], "fingerprint": fp,
                       "torsion": list(structure), "torsion_key": cells.torsion_key(structure),
                       "naive_height": naive_height(c["ainvs"]), "rank_claimed": int(c["rank_claimed"]),
                       "verification": ver,
                       "source_created_at": c.get("source_created_at"), "first_seen": c.get("first_seen")})
    if mismatches:
        raise SystemExit("!! torsion disagrees with the source for:\n  " + "\n  ".join(mismatches))
    # rank claimed descending, then ICARM before Dujella before the LMFDB, then numeric id: the first
    # match in the ledger names the known curve in KNOWN_CURVE, and a copy of an ICARM curve on a
    # Dujella page must not rename it or push the ICARM entry into duplicate_of_earlier
    curves.sort(key=lambda k: (-k["rank_claimed"], _SOURCE_ORDER.get(k["source"], 9), known_curves._id_key(k["id"])))
    # duplicate_of_earlier is defined in output order (set here only)
    seen = set()
    for k in curves:
        key = tuple(k["fingerprint"])
        k["duplicate_of_earlier"] = key in seen
        seen.add(key)
    ledger = {"version": 2, "frozen_at": raw["fetched_at"], "primes_upto": fp_mod.PRIME_BOUND, "curves": curves}
    known_curves.validate_ledger(ledger)
    return ledger


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", type=pathlib.Path, default=RAW)
    ap.add_argument("--out", type=pathlib.Path, default=OUT)
    ap.add_argument("--jobs", type=int, default=8)
    args = ap.parse_args()
    from certify_known import gp_bin
    gp = gp_bin()
    previous = json.loads(args.out.read_text()) if args.out.exists() else None
    ledger = build_ledger(load_raw(args.raw), previous, gp, args.jobs)
    args.out.write_text(dump_per_line(ledger, "curves"))
    n = len(ledger["curves"])
    print(f"wrote {args.out}: {n} curves, {sum(c['duplicate_of_earlier'] for c in ledger['curves'])} fingerprint duplicates, "
          f"{sum(c['verification']['status'] == 'certified' for c in ledger['curves'])} certified")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
