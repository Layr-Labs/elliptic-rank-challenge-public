"""The board: every submission/board/*.ec verified, deduplicated, placed in its cell, gated
against the parent board, and summed (docs/verifier-architecture.md §4, §5)."""

from __future__ import annotations

import dataclasses
import os
import pathlib
import time
from dataclasses import dataclass

from . import cells, files as files_mod, fingerprint as fp_mod, records as records_mod
from .result import ENVIRONMENT_CODES, CurveResult, ErrorCode, VerifyError, internal_error_detail
from .score import SCORING_VERSION
from .verify import VerifyConfig, verify_text

WALL_BUDGET_S = 1200.0            # whole board, parent re-verification included (docs/verifier-architecture.md §5); the harness and the CLI share it


@dataclass(frozen=True)
class CellState:
    key: str
    torsion: list
    rank: int
    known_record: dict | None
    baseline: float | None
    board_best_file: str | None
    board_best_H: float | None
    contribution: float


@dataclass(frozen=True)
class BoardResult:
    curves: list
    cells: dict
    board_total: float
    improvements: list
    parent_total: float | None
    improvement: float | None
    parent_gate: str
    ignored_parent_files: list
    records_version: int
    known_frozen_at: str
    pari_version: str

    def to_json(self) -> dict:
        curves = []
        for name, r in self.curves:
            st = self.cells[r.cell]
            d = r.to_json()
            d["file"] = name
            d["cell_record"] = st.known_record
            d["contribution"] = cells.curve_contribution(r.certified_rank, r.naive_height, st.baseline)
            d["holds_cell"] = st.board_best_file == name and st.contribution > 0.0
            curves.append(d)
        cells_json = {k: {"torsion": c.torsion, "rank": c.rank, "scored": c.rank >= cells.MIN_SCORED_RANK,
                          "known_record": c.known_record, "baseline": c.baseline,
                          "board_best": None if c.board_best_file is None else {"file": c.board_best_file, "naive_height": c.board_best_H},
                          "contribution": c.contribution}
                      for k, c in sorted(self.cells.items())}
        return {"scoring_version": SCORING_VERSION, "records_version": self.records_version,
                "known_frozen_at": self.known_frozen_at, "pari_version": self.pari_version,
                "parent_gate": self.parent_gate, "board_total": self.board_total,
                "parent_total": self.parent_total, "improvement": self.improvement,
                "curves": curves, "cells": cells_json, "improvements": self.improvements,
                "ignored_parent_files": self.ignored_parent_files}


def _verify_named(name: str, text: str, cfg: VerifyConfig) -> CurveResult:
    try:
        return verify_text(text, cfg)
    except VerifyError as e:
        raise VerifyError(e.code, f"board/{name}: {e.message}")
    except Exception as e:                     # garbage from gp (KeyError/ValueError reading its output): fail closed, named
        raise VerifyError(ErrorCode.GP_ERROR, f"board/{name}: {internal_error_detail(e)}")


def _best_per_cell(results: list) -> dict:
    """[(file, CurveResult)] -> cell -> (file, H) of the smallest H (ties: first file name)."""
    best: dict = {}
    for name, r in results:
        cur = best.get(r.cell)
        if cur is None or r.naive_height < cur[1]:
            best[r.cell] = (name, r.naive_height)
    return best


def _contributions(best: dict, table: dict) -> dict:
    """cell -> contribution of its smallest curve, (file, H) as in _best_per_cell; rank-0 cells give 0."""
    return {key: cells.curve_contribution(cells.cell_rank(key), H, records_mod.baseline(table, key))
            for key, (_, H) in best.items()}


def verify_board(files: dict, parent_files: dict | None, config: VerifyConfig | None = None,
                 records_table: dict | None = None, budget_check=None,
                 ignored_parent: list | None = None) -> BoardResult:
    """Verify, deduplicate, place and sum the board; gate it against `parent_files` if given.

    The parent board is verified first with the candidate's own settings except the conductor,
    which is off for files the candidate does not carry. A promoted file verified when its
    frontier was accepted and, with a frozen table and ledger and a pinned gp, keeps verifying;
    one that fails now (a gp failure, or a curve or point code, including KNOWN_CURVE after a
    ledger rebuild) means the environment changed under the frontier, and the run is rejected
    rather than the cell silently unheld. A parent file whose text is also a candidate file is
    verified once; if that verification fails on a per-file gate, the rejection is the
    candidate's, under the candidate's name. Every other gate still runs on every file.

    A parent cell is held when its contribution is positive. It is satisfied by a candidate curve
    there no larger than the parent's, or vacated by the parent's holder (same fingerprint: the
    curve up to isogeny and twist) being on the candidate board certified at a higher rank, in
    which case the board total must not fall. That total check is global, so one destination
    cell's gain can never be credited to two vacated cells.
    Keeping the old copy as well would be DUPLICATE_CURVE. `ignored_parent` entries (parent
    entries the harness loader skipped: bad names, non-regular files, oversize) are listed in
    `ignored_parent_files`; the same loader would have rejected the parent submission outright,
    so such an entry never held a cell."""
    cfg = config or VerifyConfig()
    parent_cfg = dataclasses.replace(cfg, conductor_enabled=False)
    tick = budget_check or (lambda: None)
    table = records_table if records_table is not None else records_mod.load()
    names = sorted(files)
    if len(names) > cells.MAX_BOARD_FILES:
        raise VerifyError(ErrorCode.TOO_MANY_CURVES, f"{len(names)} curve files; at most {cells.MAX_BOARD_FILES}")

    candidate_by_text: dict = {}
    for n in names:
        candidate_by_text.setdefault(files[n], n)
    parent_results, ignored = [], list(ignored_parent or [])
    for n in sorted(parent_files or {}):
        tick()
        text = parent_files[n]
        shared = candidate_by_text.get(text)
        try:
            parent_results.append((n, verify_text(text, cfg if shared is not None else parent_cfg)))
        except VerifyError as e:
            if e.code in ENVIRONMENT_CODES:
                raise VerifyError(e.code, f"re-verifying parent board/{n}: {e.message}")
            if shared is not None:             # the candidate's own file: its rejection, once, under its name
                raise VerifyError(e.code, f"board/{shared}: {e.message}")
            # A promoted file verified when its frontier was accepted and, with a frozen table and a
            # pinned gp, keeps verifying. Failing now means the environment changed under the frontier
            # (a ledger rebuild, a verifier change); its cell cannot be checked, so fail closed.
            raise VerifyError(ErrorCode.BOARD_REGRESSION,
                              f"previous board/{n} no longer verifies ({e.code.value}: {e.message}); the cell it held "
                              "cannot be checked and the frontier board needs a maintainer")
        except Exception as e:                 # not a verdict on the file: this run's gp output was unreadable
            raise VerifyError(ErrorCode.GP_ERROR, f"re-verifying parent board/{n}: {internal_error_detail(e)}")
    verified_by_text = {}
    for n, r in parent_results:
        verified_by_text.setdefault(parent_files[n], r)

    # Each distinct text is verified once (a parent copy, or the first candidate carrying it); two
    # candidate files with the same bytes are the same curve and DUPLICATE_CURVE below names both.
    results = []
    for n in names:
        r = verified_by_text.get(files[n])
        if r is None:
            tick()
            r = _verify_named(n, files[n], cfg)
            verified_by_text[files[n]] = r
        results.append((n, r))
    result_by_file = dict(results)
    tick()                                     # the last verification counts against the budget too
    for i in range(len(results)):
        for j in range(i + 1, len(results)):
            if fp_mod.matches(results[i][1].fingerprint, results[j][1].fingerprint):
                raise VerifyError(ErrorCode.DUPLICATE_CURVE,
                                  f"board/{results[i][0]} and board/{results[j][0]} are the same curve up to "
                                  "isogeny or quadratic twist; a curve counts once")

    best = _best_per_cell(results)
    contribs = _contributions(best, table)
    keys = set(best) | set(table["cells"])
    states: dict = {}
    for key in sorted(keys):
        rec = records_mod.record(table, key)
        base = records_mod.baseline(table, key)
        if key in best:
            fname, H = best[key]
            r = result_by_file[fname]
            torsion, rank = r.torsion, r.certified_rank
            contrib = contribs[key]
        else:
            c = table["cells"][key]
            fname, H, torsion, rank, contrib = None, None, c["torsion"], c["rank"], 0.0
        states[key] = CellState(key, list(torsion), int(rank), rec, base, fname, H, contrib)

    total = cells.board_total(contribs)
    parent_gate, parent_best, parent_contribs, parent_total = "skipped", {}, {}, None
    if parent_files is not None:
        parent_gate = "checked"
        parent_by_name = dict(parent_results)
        parent_best = _best_per_cell(parent_results)
        parent_contribs = _contributions(parent_best, table)
        parent_total = cells.board_total(parent_contribs)
        vacated = []
        for key, (pf, pH) in sorted(parent_best.items()):
            if parent_contribs[key] <= 0.0:
                continue                       # not held: rank 0, or at or above the frozen baseline
            cand = best.get(key)
            if cand is not None and cand[1] <= pH:
                continue
            # The only way to vacate a held cell: the holder itself (same fingerprint, so the curve
            # up to isogeny and twist) is on the candidate board certified at a higher rank.
            holder = parent_by_name[pf]
            if any(r.certified_rank > holder.certified_rank and fp_mod.matches(r.fingerprint, holder.fingerprint)
                   for _, r in results):
                vacated.append(key)
                continue
            have = "no curve" if cand is None else f"H = {cand[1]:.6f}"
            raise VerifyError(ErrorCode.BOARD_REGRESSION,
                              f"cell {key}: the previous board holds it with H = {pH:.6f} ({pf}); this board has {have}")
        # Moving up is allowed only while the board total does not fall. Checked once, globally: a gain
        # in one destination cell cannot be credited to two vacated cells, and whether the moved curve
        # or a new smaller one is best in the destination does not matter.
        if vacated and total < parent_total:
            raise VerifyError(ErrorCode.BOARD_REGRESSION,
                              f"cells {', '.join(vacated)} were vacated by curves moving to a higher rank, but the "
                              f"board total falls from {parent_total:.6f} to {total:.6f}")

    improvements = []
    for key, (fname, H) in sorted(best.items()):
        c, pc = contribs[key], parent_contribs.get(key, 0.0)
        if c > pc:
            improvements.append({"cell": key, "file": fname, "H": H, "contribution": c,
                                 "parent_contribution": pc, "delta": c - pc})
    improvement = None if parent_total is None else total - parent_total
    pari = results[0][1].pari_version if results else (cfg.pari_version or "")
    return BoardResult(curves=results, cells=states, board_total=total, improvements=improvements,
                       parent_total=parent_total, improvement=improvement,
                       parent_gate=parent_gate, ignored_parent_files=ignored,
                       records_version=int(table["version"]), known_frozen_at=str(table["ledger_frozen_at"]),
                       pari_version=pari)


def verify_board_dirs(board_dir: pathlib.Path, parent_dir: pathlib.Path | None, cfg: VerifyConfig,
                      records_table: dict | None = None) -> BoardResult:
    """The one wiring the harness and the CLI share, so the laptop reproduces the official verdict:
    the platform's byte cap over the editable tree (the nearest ancestor named `submission`, i.e.
    submission/ for submission/board; a board under no such directory is capped on itself), the strict
    loader for the board (a missing directory is the empty board), the lenient loader for the parent
    (a set parent path that is not a directory is a rejection, never an empty parent that passes the
    gate vacuously), and WALL_BUDGET_S as both a tick between files and a gp deadline."""
    t0 = time.monotonic()
    budget = WALL_BUDGET_S

    def budget_check():
        if time.monotonic() - t0 > budget:
            raise VerifyError(ErrorCode.RESOURCE_EXCEEDED, f"board verification exceeded {int(budget)} s")

    cfg = dataclasses.replace(cfg, deadline=t0 + budget)
    where = pathlib.Path(os.path.abspath(board_dir))            # "." inside submission/board is still submission/board
    files_mod.check_submission_tree(next((p for p in where.parents if p.name == "submission"), where))
    files = files_mod.load_board_dir(board_dir)
    parent, ignored = None, []
    if parent_dir is not None:
        if not parent_dir.is_dir():
            raise VerifyError(ErrorCode.SUBMISSION_NOT_REGULAR_FILE, f"parent board {parent_dir} is not a directory")
        parent, ignored = files_mod.load_board_dir_lenient(parent_dir)
    return verify_board(files, parent, cfg, records_table, budget_check=budget_check, ignored_parent=ignored)
