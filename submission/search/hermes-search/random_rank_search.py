#!/usr/bin/env python3
"""Small EC Rank screening helper.

This is deliberately a screening tool, not a verifier. It uses the pinned
PARI/GP binary in tools/pari/bin/gp to compute the torsion group and a
Mordell-Weil rank lower bound for random short Weierstrass models. Any hit must
still be written in the official .ec format with independent rational points
and checked by ec_verify before it can be submitted.
"""

from __future__ import annotations

import argparse
import random
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
GP = ROOT / "tools" / "pari" / "bin" / "gp"


def rank_line(ainvs: tuple[int, int, int, int, int], timeout: int = 20) -> str:
    a1, a2, a3, a4, a6 = ainvs
    expr = (
        f"E=ellinit([{a1},{a2},{a3},{a4},{a6}]);"
        "if(ellminimalmodel(E)==0,print(\"singular\"),"
        "T=elltors(E);R=ellrank(E,1);"
        "print(T[1],\" \",T[2],\" \",R[1],\" \",R[2]))"
    )
    try:
        out = subprocess.run(
            [str(GP), "-q"],
            input=expr + "\n",
            text=True,
            capture_output=True,
            timeout=timeout,
            check=False,
        )
    except subprocess.TimeoutExpired:
        return "timeout"
    return (out.stdout or out.stderr).strip().splitlines()[-1] if (out.stdout or out.stderr).strip() else "empty"


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--count", type=int, default=100)
    p.add_argument("--rank", type=int, default=1)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--bound", type=int, default=200)
    args = p.parse_args()

    if not GP.exists():
        raise SystemExit(f"missing pinned PARI/GP binary: {GP}")

    rng = random.Random(args.seed)
    found = 0
    for _ in range(args.count):
        a1 = rng.randint(-args.bound, args.bound)
        a2 = rng.randint(-args.bound, args.bound)
        a3 = rng.randint(-args.bound, args.bound)
        a4 = rng.randint(-args.bound, args.bound)
        a6 = rng.randint(-args.bound, args.bound)
        line = rank_line((a1, a2, a3, a4, a6))
        parts = line.split()
        if len(parts) < 3:
            continue
        try:
            torsion_order = int(parts[0])
            rank_lb = int(parts[2])
        except ValueError:
            continue
        if rank_lb >= args.rank:
            print(a1, a2, a3, a4, a6, "torsion", torsion_order, "rank_lb", rank_lb, flush=True)
            found += 1
    print(f"# screened={args.count} hits={found}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
