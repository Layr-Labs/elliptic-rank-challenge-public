#!/usr/bin/env python3
"""Step 2 of the record-store pipeline (docs/records.md §2): certify the published points of every raw
entry with our own verifier (novelty off), incrementally.

    python tools/certify_known.py [--raw tools/known_raw.json] [--jobs 8] [--force] [--limit N]

Only entries without a verification block, or whose source_updated_at moved since they were
verified, are run. Results are written back into the raw mirror; tools/build_known.py copies the
summary into the ledger. Only `certified` entries define records; a `failed` entry keeps its code.
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import datetime as dt
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

import ec_verify  # noqa: E402
from ec_verify import gp as gpmod  # noqa: E402
from ec_verify.result import ENVIRONMENT_CODES, VerifyError  # noqa: E402
from ec_verify.verify import VerifyConfig, verify_text  # noqa: E402
from fetch_known import RAW, load_raw, save_raw  # noqa: E402


def gp_bin() -> str:
    b = gpmod.find_gp()
    gpmod.require_pinned(b)
    return b


def submission_text(entry: dict) -> str:
    lines = ["E " + " ".join(str(a) for a in entry["ainvs"])]
    lines += [f"P {x} {y}" for x, y in entry["points"]]
    return "\n".join(lines) + "\n"


def needs_certification(entry: dict) -> bool:
    v = entry.get("verification")
    return v is None or v.get("source_updated_at") != entry.get("source_updated_at")


def certify_entry(entry: dict, gp: str, now: str, pari_version: str | None = None) -> dict:
    base = {"certified_rank": None, "code": None, "message": None, "verified_at": now,
            "pari_version": pari_version if pari_version is not None else gpmod.gp_version(gp),
            "verifier_version": ec_verify.__version__,
            "source_updated_at": entry.get("source_updated_at")}
    if not entry["points"]:
        return {"status": "unverified", **base}
    try:
        r = verify_text(submission_text(entry), VerifyConfig(gp_bin=gp, pari_version=base["pari_version"],
                                                             novelty_enabled=False, conductor_enabled=False))
    except VerifyError as e:
        if e.code in ENVIRONMENT_CODES:        # a gp timeout or crash says nothing about the curve: stop, do not record
            raise
        return {"status": "failed", **base, "code": e.code.value, "message": e.message[:300]}
    return {"status": "certified", **base, "certified_rank": r.certified_rank}


def certify_all(raw: dict, gp: str, jobs: int, now: str, force: bool = False, limit: int | None = None) -> int:
    todo = [c for c in raw["curves"] if force or needs_certification(c)]
    if limit is not None:
        todo = todo[:limit]
    pari_version = gpmod.gp_version(gp)                  # once per run, not once per curve
    with cf.ThreadPoolExecutor(max(1, jobs)) as ex:      # one gp process per curve; threads wait
        for c, v in zip(todo, ex.map(lambda c: certify_entry(c, gp, now, pari_version), todo)):
            c["verification"] = v
    return len(todo)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", type=pathlib.Path, default=RAW)
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--limit", type=int)
    args = ap.parse_args()
    raw = load_raw(args.raw)
    now = dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")
    n = certify_all(raw, gp_bin(), args.jobs, now, args.force, args.limit)
    save_raw(args.raw, raw)
    by = {}
    for c in raw["curves"]:
        by[c.get("verification", {}).get("status", "none")] = by.get(c.get("verification", {}).get("status", "none"), 0) + 1
    print(f"certified {n} this run; statuses: {by}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
