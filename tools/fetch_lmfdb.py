#!/usr/bin/env python3
"""Maintainer tool: the LMFDB's curves with their generators, from Cremona's ecdata dump.

The LMFDB web API rate-limits bulk queries behind a CAPTCHA; its elliptic-curve tables over Q are
the files of https://github.com/JohnCremona/ecdata (conductor <= 500000, every rank proven). The
`allgens` files hold one curve per line:

    N iso num [a1,a2,a3,a4,a6] rank [torsion structure] [x:y:z] ... (rank infinite-order generators, then torsion generators)

The LMFDB's small curves can only matter as records in the low-rank cells, so this tool keeps the
smallest few curves (by naive height) of every (torsion structure, rank >= 1) and writes them with
their generators to tools/lmfdb_curves.json; tools/fetch_known.py mirrors that file as the `lmfdb`
source (ids `lmfdb-<Cremona label>`, torsion published). Rank-0 curves are never records and are
not kept.

    python tools/fetch_lmfdb.py --ecdata DIR [--download] [--per-cell 3] [--out tools/lmfdb_curves.json]

--download fetches the 50 allgens files into DIR first (about 180 MiB, raw.githubusercontent.com).
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import pathlib
import subprocess
import sys
from fractions import Fraction

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from ec_verify.cells import MAZUR_GROUPS  # noqa: E402
from ec_verify.curve import naive_height  # noqa: E402

OUT = ROOT / "tools" / "lmfdb_curves.json"
RAW_BASE = "https://raw.githubusercontent.com/JohnCremona/ecdata/master/allgens/"
LIST_URL = "https://api.github.com/repos/JohnCremona/ecdata/contents/allgens"


def _rat(s: str) -> str:
    f = Fraction(s)
    return str(f.numerator) if f.denominator == 1 else f"{f.numerator}/{f.denominator}"


def parse_line(line: str) -> dict:
    """One allgens line -> label, conductor, ainvs, rank, torsion (descending, as gp and the
    verifier write it), gens (the infinite-order generators only, affine, 'num/den')."""
    f = line.split()
    if len(f) < 6:
        raise ValueError(f"short allgens line: {line[:80]}")
    n, iso, num = int(f[0]), f[1], int(f[2])
    ainvs = [int(a) for a in f[3].strip("[]").split(",")]
    if len(ainvs) != 5:
        raise ValueError(f"{n}{iso}{num}: ainvs {f[3]} is not five integers")
    rank = int(f[4])
    tors = sorted((int(x) for x in f[5].strip("[]").split(",") if x), reverse=True)
    if tuple(tors) not in MAZUR_GROUPS:
        raise ValueError(f"{n}{iso}{num}: torsion {f[5]} is not a Mazur group")
    gens = []
    for g in f[6:6 + rank]:
        x, y, z = (Fraction(int(c)) for c in g.strip("[]").split(":"))
        if z == 0:
            raise ValueError(f"{n}{iso}{num}: generator at infinity")
        gens.append([_rat(str(x / z)), _rat(str(y / z))])
    if len(gens) != rank:
        raise ValueError(f"{n}{iso}{num}: rank {rank} but {len(gens)} generators")
    return {"label": f"{n}{iso}{num}", "conductor": n, "ainvs": ainvs, "rank": rank, "torsion": tors, "gens": gens}


def select(lines, per_cell: int = 3) -> list:
    """The per_cell smallest curves (naive height, then conductor and label) of every
    (torsion structure, rank >= 1), in conductor order."""
    sys.set_int_max_str_digits(0)      # some generators in ecdata have thousands of digits
    best: dict = {}
    for line in lines:
        line = line.strip()
        if not line:
            continue
        c = parse_line(line)
        if c["rank"] < 1:
            continue
        c["naive_height"] = naive_height(c["ainvs"])
        key = (tuple(c["torsion"]), c["rank"])
        cell = best.setdefault(key, [])
        cell.append(c)
        cell.sort(key=lambda d: (d["naive_height"], d["conductor"], d["label"]))
        del cell[per_cell:]
    out = [c for cell in best.values() for c in cell]
    out.sort(key=lambda d: (d["conductor"], d["label"]))
    return out


def download(into: pathlib.Path) -> None:
    into.mkdir(parents=True, exist_ok=True)
    listing = subprocess.run(["curl", "-fsSL", "-m", "60", LIST_URL], capture_output=True, text=True, check=True).stdout
    names = sorted(e["name"] for e in json.loads(listing) if e["name"].startswith("allgens."))
    for name in names:
        dest = into / name
        if dest.exists() and dest.stat().st_size > 0:
            continue
        subprocess.run(["curl", "-fsSL", "-m", "600", RAW_BASE + name, "-o", str(dest)], check=True)
        print(f"fetched {name}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ecdata", type=pathlib.Path, required=True, help="directory holding the allgens.* files")
    ap.add_argument("--download", action="store_true", help="fetch the allgens files into --ecdata first")
    ap.add_argument("--per-cell", type=int, default=3)
    ap.add_argument("--out", type=pathlib.Path, default=OUT)
    args = ap.parse_args()
    if args.download:
        download(args.ecdata)
    files = sorted(args.ecdata.glob("allgens.*"))
    if not files:
        raise SystemExit(f"!! no allgens.* files under {args.ecdata}")

    def lines():
        for f in files:
            with f.open() as fh:
                yield from fh
    chosen = select(lines(), args.per_cell)
    doc = {"source": "JohnCremona/ecdata allgens (LMFDB elliptic curves over Q, conductor <= 500000)",
           "fetched_at": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
           "files": len(files), "per_cell": args.per_cell, "curves": chosen}
    args.out.write_text(json.dumps(doc, indent=1) + "\n")
    cells = {(tuple(c["torsion"]), c["rank"]) for c in chosen}
    print(f"wrote {args.out}: {len(chosen)} curves over {len(cells)} (torsion, rank) cells from {len(files)} files")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
