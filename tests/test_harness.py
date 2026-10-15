"""Harness e2e through the REAL entrypoint (subprocess) over the board directory. Every rejection
asserts nonzero exit AND no score file."""

import json
import os
import pathlib
import subprocess
import sys

import pytest

from conftest import needs_gp
from util import ROOT, base5, base5_dependent, gens_5077a

pytestmark = needs_gp


def _run(tmp_path, files=None, parent=None, board_entry=None, timeout=1800, parent_fifo=None, extra_env=None):
    """files: {name: text} written under the repo's submission board dir. parent: {name: text}
    exported to a parent dir and passed via EC_PARENT_BOARD_DIR. board_entry: ('symlink', target) |
    ('dir',) | ('fifo',) to plant a non-regular entry named weird.ec."""
    repo = tmp_path / "repo"
    board = repo / "submission" / "board"
    board.mkdir(parents=True)
    for name, text in (files or {}).items():
        (board / name).write_text(text, encoding="ascii")
    if board_entry:
        target = board / "weird.ec"
        if board_entry[0] == "symlink":
            os.symlink(board_entry[1], target)
        elif board_entry[0] == "dir":
            target.mkdir()
        elif board_entry[0] == "fifo":
            os.mkfifo(target)
    env = {"PYTHONHASHSEED": "0", "PYTHONPATH": str(ROOT), "PATH": os.environ.get("PATH", ""),
           "EC_GP_BIN": os.environ.get("EC_GP_BIN", "")}
    if parent is not None:
        pdir = tmp_path / "parent"
        pdir.mkdir()
        for name, text in parent.items():
            (pdir / name).write_text(text, encoding="ascii")
        if parent_fifo:
            os.mkfifo(pdir / parent_fifo)
        env["EC_PARENT_BOARD_DIR"] = str(pdir)
    env = {k: v for k, v in env.items() if v}
    env.update(extra_env or {})
    proc = subprocess.run([sys.executable, "-P", "-m", "harness.verify"], cwd=repo, env=env,
                          capture_output=True, text=True, timeout=timeout)
    score_path = repo / "score.json"
    score = json.loads(score_path.read_text()) if score_path.exists() else None
    return proc, score


def test_dependent_rejected(tmp_path):
    proc, score = _run(tmp_path, {"a.ec": base5_dependent()})
    assert proc.returncode != 0 and score is None
    assert "REJECTED: POINTS_DEPENDENT: board/a.ec:" in proc.stdout


def test_garbage_rejected(tmp_path):
    proc, score = _run(tmp_path, {"a.ec": "not a curve\n"})
    assert proc.returncode != 0 and score is None and "REJECTED: PARSE_SYNTAX" in proc.stdout


def test_missing_board_dir_is_an_empty_board(tmp_path):
    repo = tmp_path / "repo"; (repo / "submission").mkdir(parents=True)
    env = {k: v for k, v in {"PYTHONHASHSEED": "0", "PYTHONPATH": str(ROOT), "PATH": os.environ.get("PATH", ""),
                             "EC_GP_BIN": os.environ.get("EC_GP_BIN", "")}.items() if v}
    proc = subprocess.run([sys.executable, "-P", "-m", "harness.verify"], cwd=repo, env=env, capture_output=True, text=True)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    score = json.loads((repo / "score.json").read_text())
    assert score["score"] == 0.0 and score["metrics"]["curves"] == [] and score["metrics"]["parent_gate"] == "skipped"


def test_symlink_entry_rejected(tmp_path):
    real = tmp_path / "real.ec"; real.write_text(gens_5077a(), encoding="ascii")
    proc, score = _run(tmp_path, {}, board_entry=("symlink", real))
    assert proc.returncode != 0 and score is None and "SUBMISSION_NOT_REGULAR_FILE" in proc.stdout


def test_subdirectory_named_like_a_curve_is_ignored(tmp_path):
    proc, score = _run(tmp_path, {}, board_entry=("dir",))
    assert proc.returncode == 0 and score["metrics"]["curves"] == []


@pytest.mark.skipif(not hasattr(os, "mkfifo"), reason="no FIFOs on this platform")
def test_fifo_entry_rejected_without_hanging(tmp_path):
    proc, score = _run(tmp_path, {}, board_entry=("fifo",), timeout=10)
    assert proc.returncode != 0 and score is None and "SUBMISSION_NOT_REGULAR_FILE" in proc.stdout


def test_bad_file_name_rejected(tmp_path):
    proc, score = _run(tmp_path, {"my curve.ec": base5()})
    assert proc.returncode != 0 and score is None and "PARSE_SYNTAX" in proc.stdout and "Traceback" not in proc.stderr


def test_non_ec_files_are_ignored(tmp_path):
    proc, score = _run(tmp_path, {"notes.txt": "anything\n", "search.gp": "print(1)\n"})
    assert proc.returncode == 0 and score["metrics"]["curves"] == []


def test_oversize_file_rejected(tmp_path):
    proc, score = _run(tmp_path, {"a.ec": gens_5077a() + "#" + "x" * (1024 * 1024) + "\n"})
    assert proc.returncode != 0 and score is None and "RESOURCE_EXCEEDED" in proc.stdout


def test_large_legal_file_passes_the_size_gate(tmp_path):
    # 600 KB of comments plus a dependent pair: must reach the certificate, not the size cap.
    dep = base5_dependent().splitlines()
    text = dep[1] + "\n" + ("# " + "x" * 9000 + "\n") * 70 + "\n".join(dep[2:]) + "\n"
    assert len(text) > 512 * 1024
    proc, score = _run(tmp_path, {"a.ec": text})
    assert proc.returncode != 0 and score is None
    assert "POINTS_DEPENDENT" in proc.stdout


def test_known_curve_rejected(tmp_path):
    proc, score = _run(tmp_path, {"a.ec": gens_5077a()})
    assert proc.returncode != 0 and score is None and "REJECTED: KNOWN_CURVE: board/a.ec:" in proc.stdout


def test_long_claim_rejected_cleanly(tmp_path):
    proc, score = _run(tmp_path, {"a.ec": gens_5077a() + "#CLAIM rank=" + "9" * 5000 + "\n"})
    assert proc.returncode != 0 and score is None and "REJECTED: PARSE_SYNTAX" in proc.stdout and "Traceback" not in proc.stderr


def test_board_scores_and_is_deterministic(tmp_path):
    proc, score = _run(tmp_path, {"b5.ec": base5()})
    assert proc.returncode == 0, proc.stdout + proc.stderr
    m = score["metrics"]
    assert m["curves"][0]["certified_rank"] == 5 and m["curves"][0]["cell"] == "1/5"
    assert m["scoring_version"] == 6 and isinstance(m["records_version"], int)
    assert score["score"] >= 0.0                       # base5 is above the real rank-5 record: 0 unless the table says otherwise
    proc2, score2 = _run(tmp_path / "again", {"b5.ec": base5()})
    def core(d):
        d = json.loads(json.dumps(d))
        for c in d["metrics"]["curves"]:
            c.pop("conductor")
        return json.dumps(d, sort_keys=True)
    assert core(score) == core(score2) and score["score"] == score2["score"]


def test_parent_dir_enables_the_gate(tmp_path):
    proc, score = _run(tmp_path, {"b5.ec": base5()}, parent={})
    assert proc.returncode == 0 and score["metrics"]["parent_gate"] == "checked"
    assert score["metrics"]["parent_total"] == 0.0 and score["metrics"]["improvement"] == score["score"]
    from util import fixture
    from ec_verify import records
    proc, score = _run(tmp_path / "reg", {"b5.ec": base5()}, parent={"r6.ec": fixture("rank6_blind.ec")})
    # rank6_blind (cell 1/6, H 4.288426) holds its cell only if it undercuts that cell's baseline in the
    # live table (or the cell has none): then dropping it is a regression, otherwise it is allowed.
    base = records.baseline(records.load(), "1/6")
    if base is None or 4.288426 < base:
        assert proc.returncode != 0 and score is None
        assert "REJECTED: BOARD_REGRESSION: cell 1/6" in proc.stdout and "(r6.ec)" in proc.stdout
    else:
        assert proc.returncode == 0, proc.stdout


def test_unloadable_parent_entries_are_ignored_not_fatal(tmp_path):
    # A parent entry that cannot be loaded (FIFO, bad name) is the previous frontier's problem,
    # not this candidate's: it is skipped and listed in ignored_parent_files.
    proc, score = _run(tmp_path, {"b5.ec": base5()}, parent={"b5.ec": base5(), "bad name.ec": base5()},
                       parent_fifo="weird.ec", timeout=120)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    ign = {e["file"]: e for e in score["metrics"]["ignored_parent_files"]}
    assert set(ign) == {"weird.ec", "bad name.ec"}
    assert ign["weird.ec"]["code"] == "SUBMISSION_NOT_REGULAR_FILE" and ign["weird.ec"]["message"]
    assert ign["bad name.ec"]["code"] == "PARSE_SYNTAX"
    assert score["metrics"]["parent_gate"] == "checked"


def test_missing_parent_dir_is_rejected_not_checked_vacuously(tmp_path):
    # A set EC_PARENT_BOARD_DIR that does not exist must not load as an empty parent board: the gate
    # would then pass every deletion while score.json still reports parent_gate "checked".
    proc, score = _run(tmp_path, {"b5.ec": base5()}, extra_env={"EC_PARENT_BOARD_DIR": str(tmp_path / "nope")})
    assert proc.returncode != 0 and score is None
    assert proc.stdout.startswith("REJECTED: SUBMISSION_NOT_REGULAR_FILE:") and "Traceback" not in proc.stderr, proc.stdout


def test_gp_version_mismatch_and_missing_gp_report_their_own_codes(tmp_path):
    # The documented environment codes (docs/submitting.md §6), not "GP_ERROR: internal error: ...".
    fake = tmp_path / "gp"
    fake.write_text("#!/bin/sh\necho 2.15.0\n", encoding="ascii")
    fake.chmod(0o755)
    proc, score = _run(tmp_path, {"b5.ec": base5()}, extra_env={"EC_GP_BIN": str(fake)})
    assert proc.returncode != 0 and score is None and proc.stdout.startswith("REJECTED: GP_VERSION_MISMATCH:"), proc.stdout
    proc, score = _run(tmp_path / "b", {"b5.ec": base5()}, extra_env={"EC_GP_BIN": str(tmp_path / "missing-gp")})
    assert proc.returncode != 0 and score is None and proc.stdout.startswith("REJECTED: GP_UNAVAILABLE:"), proc.stdout


def test_required_parent_gate_fails_closed_when_no_parent_was_exported(tmp_path):
    # benchmark.yml sets EC_REQUIRE_PARENT with EC_PARENT_REF on every submissions/* ref. If the parent
    # export did not happen (no EC_PARENT_BOARD_DIR), the run must be rejected, never scored with
    # parent_gate "skipped".
    proc, score = _run(tmp_path, {"b5.ec": base5()}, extra_env={"EC_REQUIRE_PARENT": "1"})
    assert proc.returncode != 0 and score is None
    assert proc.stdout.startswith("REJECTED: BOARD_REGRESSION:") and "EC_PARENT_BOARD_DIR" in proc.stdout, proc.stdout
    proc, score = _run(tmp_path / "ok", {"b5.ec": base5()}, parent={}, extra_env={"EC_REQUIRE_PARENT": "1"})
    assert proc.returncode == 0 and score["metrics"]["parent_gate"] == "checked"


def test_submission_tree_over_the_platform_cap_is_rejected(tmp_path):
    # benchmark.json maxSubmissionBytes covers everything under submission/, not only board/*.ec.
    repo = tmp_path / "repo"
    notes = repo / "submission" / "notes"
    notes.mkdir(parents=True)
    (notes / "search.log").write_bytes(b"x" * (16 * 1024 * 1024 + 1))
    proc, score = _run(tmp_path, {"b5.ec": base5()})
    assert proc.returncode != 0 and score is None
    assert proc.stdout.startswith("REJECTED: RESOURCE_EXCEEDED:") and "submission" in proc.stdout, proc.stdout
