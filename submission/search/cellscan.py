"""Frozen cell table beside the board's holders: baseline H, holder H, contribution.

Usage: EC_GP_BIN=... python3 submission/search/cellscan.py [board_dir] [torsion_key ...]
"""
import pathlib
import sys

from ec_verify import gp as gpmod
from ec_verify import records as records_mod
from ec_verify.board import verify_board_dirs
from ec_verify.verify import VerifyConfig

ROOT = pathlib.Path(__file__).resolve().parents[2]


def main() -> None:
    args = sys.argv[1:]
    board = pathlib.Path(args.pop(0)) if args and "/" in args[0] else ROOT / "submission" / "board"
    only = set(args)
    gp_bin = gpmod.find_gp()
    cfg = VerifyConfig(gp_bin=gp_bin, pari_version=gpmod.gp_version(gp_bin), conductor_enabled=False,
                       novelty_enabled=False, require_pinned_gp=False)
    held = verify_board_dirs(board, None, cfg).to_json()["cells"]
    table = records_mod.load()
    for key, c in sorted(table["cells"].items(), key=lambda kv: (kv[1]["torsion_key"], kv[1]["rank"])):
        if only and c["torsion_key"] not in only:
            continue
        rec = c["record"]
        h = held.get(key)
        best = (h or {}).get("board_best")
        mark = f"{best['naive_height']:.6f} {h['contribution']:.3f} {best['file']}" if best else "-"
        print(f"{key:8} base={rec['naive_height']:<10} {rec['id']:<12} known={c['known_count']:<5} {mark}")


if __name__ == "__main__":
    main()
