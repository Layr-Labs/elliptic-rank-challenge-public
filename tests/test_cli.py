"""python -m ec_verify loads a board directory exactly as the harness does: the laptop must
reproduce the official run, including its rejections of bad names, non-regular files and
oversize files (docs/verifier-architecture.md §5)."""

import json
import os

import pytest

from conftest import needs_gp
from util import base5, fixture

from ec_verify.files import load_board_dir

from ec_verify.__main__ import main

pytestmark = needs_gp


def _run(capsys, path):
    rc = main([str(path), "--no-conductor", "--no-novelty"])
    return rc, capsys.readouterr().out


def test_cli_rejects_a_board_file_name_the_harness_rejects(tmp_path, capsys):
    (tmp_path / "bad name.ec").write_text(base5(), encoding="ascii")
    rc, out = _run(capsys, tmp_path)
    assert rc == 1 and out.startswith("REJECTED: PARSE_SYNTAX:"), out


def test_cli_rejects_a_symlinked_board_file(tmp_path, capsys):
    (tmp_path / "a.ec").write_text(base5(), encoding="ascii")
    os.symlink("a.ec", tmp_path / "link.ec")
    rc, out = _run(capsys, tmp_path)
    assert rc == 1 and out.startswith("REJECTED: SUBMISSION_NOT_REGULAR_FILE:"), out


def test_cli_rejects_an_oversize_board_file(tmp_path, capsys):
    (tmp_path / "big.ec").write_text(base5() + ("# " + "x" * 5000 + "\n") * 220, encoding="ascii")   # > 1 MiB, every line legal
    rc, out = _run(capsys, tmp_path)
    assert rc == 1 and out.startswith("REJECTED: RESOURCE_EXCEEDED:"), out


def test_cli_accepts_a_clean_board(tmp_path, capsys):
    (tmp_path / "b5.ec").write_text(base5(), encoding="ascii")
    rc, out = _run(capsys, tmp_path)
    assert rc == 0 and '"score"' in out


def test_cli_rejects_a_missing_parent_directory(tmp_path, capsys):
    # A mistyped --parent must not pass the gate vacuously against an empty parent board.
    (tmp_path / "b5.ec").write_text(base5(), encoding="ascii")
    with pytest.raises(SystemExit) as ei:
        main([str(tmp_path), "--parent", str(tmp_path / "nope"), "--no-conductor", "--no-novelty"])
    assert ei.value.code == 2


def test_cli_lists_unloadable_parent_entries_like_the_harness(tmp_path, capsys):
    board, parent = tmp_path / "board", tmp_path / "parent"
    board.mkdir(), parent.mkdir()
    (board / "b5.ec").write_text(base5(), encoding="ascii")
    (parent / "bad name.ec").write_text(base5(), encoding="ascii")
    rc = main([str(board), "--parent", str(parent), "--no-conductor", "--no-novelty"])
    out = json.loads(capsys.readouterr().out)
    assert rc == 0 and out["metrics"]["parent_gate"] == "checked"
    assert [e["code"] for e in out["metrics"]["ignored_parent_files"]] == ["PARSE_SYNTAX"]


def test_cli_single_file_uses_the_strict_loader(tmp_path, capsys):
    big = tmp_path / "big.ec"
    big.write_text(base5() + ("# " + "x" * 5000 + "\n") * 220, encoding="ascii")
    rc, out = _run(capsys, big)
    assert rc == 1 and out.startswith("REJECTED: RESOURCE_EXCEEDED:"), out
    with pytest.raises(SystemExit) as ei:                  # a path that does not exist is a usage error, not a verdict
        _run(capsys, tmp_path / "missing.ec")
    assert ei.value.code == 2


def test_cli_board_mode_checks_the_gp_pin_once_not_per_curve(tmp_path, monkeypatch):
    # The harness passes pari_version into the config; the CLI must do the same or a 400-file
    # board plus parent spends ~800 gp launches on `gp --version-short`.
    import ec_verify.gp as gpm
    launches = []
    real = gpm.gp_version                      # the one function that runs `gp --version-short` (require_pinned calls it)
    monkeypatch.setattr(gpm, "gp_version", lambda b: launches.append(b) or real(b))
    board, parent = tmp_path / "board", tmp_path / "parent"
    board.mkdir(), parent.mkdir()
    (board / "b5.ec").write_text(base5(), encoding="ascii")
    (parent / "a.ec").write_text(base5(), encoding="ascii")
    (parent / "r6.ec").write_text(fixture("rank6_blind.ec"), encoding="ascii")
    rc = main([str(board), "--parent", str(parent), "--no-conductor", "--no-novelty"])
    assert rc == 0 and len(launches) == 1, launches


def test_cli_refuses_parent_with_a_single_file(tmp_path, capsys):
    # --parent only means something for a board; silently ignoring it would print a result with
    # no regression gate while the user believes one ran.
    (tmp_path / "b5.ec").write_text(base5(), encoding="ascii")
    (tmp_path / "parent").mkdir()
    with pytest.raises(SystemExit) as ei:
        main([str(tmp_path / "b5.ec"), "--parent", str(tmp_path / "parent"), "--no-conductor", "--no-novelty"])
    assert ei.value.code == 2


def test_dangling_symlink_as_the_board_is_rejected_not_an_empty_board(tmp_path):
    # Path.exists() follows symlinks: a dangling one must not load as an empty board (the harness
    # would exit 0 with score 0); the symlink rejection applies before the missing-directory case.
    os.symlink(tmp_path / "nonexistent", tmp_path / "board")
    from ec_verify import ErrorCode, VerifyError
    with pytest.raises(VerifyError) as ei:
        load_board_dir(tmp_path / "board")
    assert ei.value.code is ErrorCode.SUBMISSION_NOT_REGULAR_FILE
    assert load_board_dir(tmp_path / "missing") == {}            # a plain missing directory stays an empty board


def test_cli_single_file_reports_the_cell_record_and_contribution(tmp_path, capsys):
    # One curve: what it alone would earn in its cell under the frozen table (docs/verifier-architecture.md §4), no percentages.
    f = tmp_path / "b5.ec"
    f.write_text(base5(), encoding="ascii")
    rc, out = _run(capsys, f)
    assert rc == 0, out
    d = json.loads(out)
    from ec_verify import cells, records
    t = records.load()
    assert d["cell"] == "1/5" and d["cell_record"] == records.record(t, "1/5")
    assert d["contribution"] == pytest.approx(cells.contribution(d["naive_height"], records.baseline(t, "1/5")))
    assert "pct_vs_record" not in d and "contribution_vs_baseline" not in d


def test_cli_board_mode_enforces_the_same_budget_as_the_harness(tmp_path, capsys, monkeypatch):
    # "Your laptop reproduces the score exactly": a board the official run rejects for exceeding the
    # 1200 s budget must be rejected locally too, not scored.
    import ec_verify.board as bmod
    monkeypatch.setattr(bmod, "WALL_BUDGET_S", -1.0)
    (tmp_path / "b5.ec").write_text(base5(), encoding="ascii")
    rc, out = _run(capsys, tmp_path)
    assert rc == 1 and out.startswith("REJECTED: RESOURCE_EXCEEDED: board verification exceeded"), out


def test_cli_rejects_a_board_over_the_platform_byte_cap(tmp_path, capsys, monkeypatch):
    # Yukon refuses a submission above benchmark.json's maxSubmissionBytes before the harness runs;
    # the same total is enforced here so the laptop sees a REJECTED code instead.
    import ec_verify.files as fmod
    monkeypatch.setattr(fmod, "MAX_BOARD_BYTES", len(base5()) + 10)
    (tmp_path / "a.ec").write_text(base5(), encoding="ascii")
    (tmp_path / "b.ec").write_text(fixture("rank6_blind.ec"), encoding="ascii")
    rc, out = _run(capsys, tmp_path)
    assert rc == 1 and out.startswith("REJECTED: RESOURCE_EXCEEDED:") and "totals" in out, out


def test_cli_refuses_a_nonexistent_path_instead_of_scoring_an_empty_board(tmp_path):
    # The harness scores a missing board directory as the empty board (a submission with no files is
    # legal); the CLI must not turn a typo into a well-formed score of 0 with exit 0 that a search loop
    # takes for a verdict on its board. Exit 2 (usage error), distinct from a REJECTED verdict (exit 1).
    with pytest.raises(SystemExit) as ei:
        main([str(tmp_path / "submisson" / "board"), "--no-conductor", "--no-novelty"])
    assert ei.value.code == 2


def test_cli_loader_rejections_carry_the_board_file_prefix(tmp_path, capsys):
    # CLAUDE.md: "the detail starts with board/<file>:" — for loader-level failures too.
    (tmp_path / "big.ec").write_bytes(b"#" + b"x" * (1024 * 1024))
    rc, out = _run(capsys, tmp_path)
    assert rc == 1 and out.startswith("REJECTED: RESOURCE_EXCEEDED: board/big.ec: "), out
    (tmp_path / "big.ec").unlink()
    (tmp_path / "bad name.ec").write_text(base5(), encoding="ascii")
    rc, out = _run(capsys, tmp_path)
    assert rc == 1 and out.startswith("REJECTED: PARSE_SYNTAX: board/bad name.ec: "), out


def test_cli_rejects_a_dotfile_board_name_like_the_workflow_surface_check(tmp_path, capsys):
    # benchmark.yml rejects every dotfile under submission/ before the harness runs; the loader
    # must agree so the laptop never scores a board the official run refuses without a code.
    (tmp_path / ".wip.ec").write_text(base5(), encoding="ascii")
    rc, out = _run(capsys, tmp_path)
    assert rc == 1 and out.startswith("REJECTED: PARSE_SYNTAX: board/.wip.ec: "), out


def test_cli_reports_an_internal_error_as_a_rejection_not_a_traceback(tmp_path, capsys, monkeypatch):
    from ec_verify import __main__ as cli
    (tmp_path / "b5.ec").write_text(base5(), encoding="ascii")

    def boom(*a, **k):
        raise KeyError("AGREE")
    monkeypatch.setattr(cli.board_mod, "verify_board_dirs", boom)
    rc, out = _run(capsys, tmp_path)
    assert rc == 1 and out.startswith("REJECTED: GP_ERROR: internal error: KeyError"), out
