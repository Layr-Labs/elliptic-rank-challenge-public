"""Yukon evaluator: grade every submission/board/*.ec, write score.json.

Run as `python -E -s -m harness.verify` from the repo root (benchmark.sh runs it
inside a sandbox with EC_SCORE_DIR pointing at scratch and EC_GP_BIN at the
pinned gp). The board is data only: our parser and the pinned PARI/GP are the
only things that ever read it; nothing in it is executed. EC_PARENT_BOARD_DIR,
when set, names the previous board directory and enables the regression gate;
EC_REQUIRE_PARENT, when set, makes a run without it a rejection (fail closed).

score.json is written only on full success; every rejection prints
`REJECTED: <CODE>: <detail>`, exits 1, and writes no score.
"""

from __future__ import annotations

import json
import math
import os
import pathlib
import sys

from ec_verify import gp as gpmod
from ec_verify import board as board_mod
from ec_verify.result import VerifyError, internal_error_detail
from ec_verify.score import yukon_score
from ec_verify.verify import VerifyConfig, gp_guard

ROOT = pathlib.Path.cwd()
BOARD_DIR = ROOT / "submission" / "board"
SCORE_PATH = pathlib.Path(os.environ.get("EC_SCORE_DIR", str(ROOT))) / "score.json"


class Reject(Exception):
    pass


def rejecting(fn, *args, **kwargs):
    """Run an ec_verify step, surfacing its VerifyError as Reject (the harness's one interface)."""
    try:
        return fn(*args, **kwargs)
    except VerifyError as e:
        raise Reject(f"{e.code.value}: {e.message}")


def main() -> None:
    try:
        SCORE_PATH.unlink()
    except OSError:
        pass
    pdir = os.environ.get("EC_PARENT_BOARD_DIR")            # benchmark.sh always exports the directory when set
    if os.environ.get("EC_REQUIRE_PARENT") and not pdir:
        # benchmark.yml requires the gate on every submissions/* ref (EC_PARENT_REF and EC_REQUIRE_PARENT
        # set together; baseline and maintainer runs set neither). No parent here means the plumbing broke:
        # fail closed rather than score a board that nobody gated against the frontier. This is the single
        # enforcement point.
        raise Reject("BOARD_REGRESSION: the regression gate is required for this run (EC_REQUIRE_PARENT is set) but no "
                     "parent board was exported (EC_PARENT_BOARD_DIR unset); needs a maintainer (benchmark.yml, "
                     "benchmark.sh, EC_PARENT_REF)")
    gp_bin = os.environ.get("EC_GP_BIN") or rejecting(gp_guard, gpmod.find_gp)
    pari_version = rejecting(gp_guard, gpmod.require_pinned, gp_bin)       # once per run; verify_text reuses it
    cfg = VerifyConfig(gp_bin=gp_bin, pari_version=pari_version)
    # Loaders, budget tick and gp deadline are ec_verify.board.verify_board_dirs: the CLI runs the same.
    board = rejecting(board_mod.verify_board_dirs, BOARD_DIR, pathlib.Path(pdir) if pdir else None, cfg)
    score = yukon_score(board)
    if not math.isfinite(score):
        raise Reject("GP_ERROR: computed score is not finite")
    payload = {"score": score, "metrics": board.to_json()}
    SCORE_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(SCORE_PATH, "w", encoding="ascii", newline="\n") as fh:
        json.dump(payload, fh, sort_keys=True)
    held = sum(1 for c in board.cells.values() if c.contribution > 0.0)
    print(f"score {score} (curves={len(board.curves)}, cells_held={held}, improvement={board.improvement}, "
          f"records_version={board.records_version})", flush=True)


if __name__ == "__main__":
    try:
        main()
    except Reject as e:
        print(f"REJECTED: {e}", flush=True)
        sys.exit(1)
    except Exception as e:  # never a traceback as the interface; still exit 1, no score
        print(f"REJECTED: GP_ERROR: {internal_error_detail(e)}", flush=True)
        sys.exit(1)
