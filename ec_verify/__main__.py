"""Local CLI.

    python -m ec_verify submission/board            # the board: {"score", "metrics"} as the harness writes it
    python -m ec_verify path/to/curve.ec            # one curve: its result, cell record and contribution
    python -m ec_verify submission/board --parent /path/to/parent/board   # enable the regression gate
"""

from __future__ import annotations

import argparse
import json
import pathlib
import sys

from . import board as board_mod, cells, gp as gpmod, records
from .files import strict_load_text
from .result import VerifyError, internal_error_detail
from .score import yukon_score
from .verify import VerifyConfig, gp_guard, verify_text


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="python -m ec_verify")
    ap.add_argument("path", help="a .ec file or a board directory")
    ap.add_argument("--parent", help="previous board directory (enables BOARD_REGRESSION and the improvement against it)")
    ap.add_argument("--no-conductor", action="store_true")
    ap.add_argument("--no-novelty", action="store_true", help="skip the known-curve gate (search loops only)")
    ap.add_argument("--allow-gp-version-mismatch", action="store_true")
    args = ap.parse_args(argv)
    p = pathlib.Path(args.path)
    if not p.exists() and not p.is_symlink():
        # The harness treats a missing submission/board as the empty board; a CLI typo must not become a
        # well-formed score of 0 with exit 0 that a search loop takes for a verdict on its board.
        ap.error(f"{args.path}: no such file or directory")
    board_mode = p.is_dir() or p.suffix != ".ec"
    parent_dir = pathlib.Path(args.parent) if args.parent else None
    if parent_dir is not None and not parent_dir.is_dir():
        ap.error(f"--parent {args.parent}: not a directory (a missing parent would pass the gate vacuously)")
    if parent_dir is not None and not board_mode:
        ap.error(f"--parent applies to a board directory; {args.path} is not one (a single .ec file has no regression gate)")
    try:
        # The pin is checked once here, as the harness does; verify_text then launches no gp for it per curve.
        gp_bin = gp_guard(gpmod.find_gp)
        pari_version = gp_guard(gpmod.gp_version if args.allow_gp_version_mismatch else gpmod.require_pinned, gp_bin)
        cfg = VerifyConfig(gp_bin=gp_bin, pari_version=pari_version, conductor_enabled=not args.no_conductor,
                           novelty_enabled=not args.no_novelty, require_pinned_gp=not args.allow_gp_version_mismatch)
        if board_mode:
            # Exactly the harness's wiring: tree cap, loaders, budget tick, gp deadline (ec_verify.board.verify_board_dirs).
            b = board_mod.verify_board_dirs(p, parent_dir, cfg)
            out = {"score": yukon_score(b), "metrics": b.to_json()}
        else:
            r = verify_text(strict_load_text(p), cfg)
            t = records.load()
            rec = records.record(t, r.cell)
            out = r.to_json()
            out["cell_record"] = rec
            out["contribution"] = cells.curve_contribution(r.certified_rank, r.naive_height, records.baseline(t, r.cell))
    except VerifyError as e:
        print(f"REJECTED: {e.code.value}: {e.message}")
        return 1
    except Exception as e:                     # the harness's interface, never a traceback
        print(f"REJECTED: GP_ERROR: {internal_error_detail(e)}")
        return 1
    print(json.dumps(out, indent=1, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
