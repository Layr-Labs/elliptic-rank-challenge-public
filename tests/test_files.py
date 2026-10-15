"""The platform's submission cap covers the whole editable tree, so the local check must too."""

import os
import pathlib

import pytest

from ec_verify import ErrorCode, VerifyError
from ec_verify.files import submission_tree_bytes, check_submission_tree


def test_tree_bytes_counts_regular_files_recursively_without_following_symlinks(tmp_path):
    (tmp_path / "board").mkdir()
    (tmp_path / "board" / "a.ec").write_bytes(b"x" * 10)
    (tmp_path / "notes").mkdir()
    (tmp_path / "notes" / "log.txt").write_bytes(b"y" * 20)
    os.symlink(tmp_path / "notes" / "log.txt", tmp_path / "link.txt")      # a symlink is not counted
    assert submission_tree_bytes(tmp_path) == 30


def test_check_rejects_a_tree_over_the_cap_and_names_the_bytes(tmp_path):
    (tmp_path / "big.bin").write_bytes(b"z" * 101)
    check_submission_tree(tmp_path, max_bytes=101)
    with pytest.raises(VerifyError) as ei:
        check_submission_tree(tmp_path, max_bytes=100)
    assert ei.value.code is ErrorCode.RESOURCE_EXCEEDED and "101" in ei.value.message


def test_check_of_a_missing_tree_is_zero_bytes(tmp_path):
    assert submission_tree_bytes(tmp_path / "nope") == 0


def test_tree_cap_counts_utf8_bytes_not_code_points(tmp_path):
    # 'é' is one code point and two bytes: "E 1é" is 4 code points and 5 bytes. A cap of 4 bytes
    # must reject it (counting code points would let it through).
    (tmp_path / "a.ec").write_text("E 1é", encoding="utf-8")
    with pytest.raises(VerifyError) as ei:
        check_submission_tree(tmp_path, max_bytes=4)
    assert ei.value.code is ErrorCode.RESOURCE_EXCEEDED and "5 bytes" in ei.value.message


def test_board_dirs_cap_submission_for_submission_board_else_the_directory(tmp_path, monkeypatch):
    # a board directly under submission/ -> the cap covers submission/ (what the platform measures); a board under no
    # directory named submission is capped on itself. One rule, in the wiring the harness and the CLI share.
    from ec_verify import files
    from ec_verify.board import verify_board_dirs
    from ec_verify.verify import VerifyConfig
    sub = tmp_path / "submission"
    (sub / "board").mkdir(parents=True)
    (sub / "notes.txt").write_bytes(b"n" * 50)
    monkeypatch.setattr(files, "MAX_BOARD_BYTES", 40)
    with pytest.raises(VerifyError) as ei:
        verify_board_dirs(sub / "board", None, VerifyConfig())
    assert ei.value.code is ErrorCode.RESOURCE_EXCEEDED and "submission/ totals 50 bytes" in ei.value.message
    other = tmp_path / "elsewhere"
    other.mkdir()
    (other / "notes.txt").write_bytes(b"n" * 50)
    with pytest.raises(VerifyError) as ei:
        verify_board_dirs(other, None, VerifyConfig())
    assert "elsewhere/ totals 50 bytes" in ei.value.message


def test_board_dirs_cap_root_does_not_depend_on_how_the_path_is_spelled(tmp_path, monkeypatch):
    # `python -m ec_verify .` from inside the board directory must cap the same tree as the harness does
    # from the repo root: the rule looks at the normalised absolute path, not the literal argument.
    from ec_verify import files
    from ec_verify.board import verify_board_dirs
    from ec_verify.verify import VerifyConfig
    sub = tmp_path / "submission"
    (sub / "board").mkdir(parents=True)
    (sub / "notes.txt").write_bytes(b"n" * 50)
    monkeypatch.setattr(files, "MAX_BOARD_BYTES", 40)
    monkeypatch.chdir(sub / "board")
    with pytest.raises(VerifyError) as ei:
        verify_board_dirs(pathlib.Path("."), None, VerifyConfig())
    assert "submission/ totals 50 bytes" in ei.value.message


def test_tree_bytes_fails_closed_on_an_unreadable_subdirectory(tmp_path):
    # An unreadable directory under submission/ hides bytes the platform will count: refuse to measure
    # rather than report an undercount that the platform then rejects.
    if os.geteuid() == 0:
        pytest.skip("root reads everything")
    hidden = tmp_path / "notes"
    hidden.mkdir()
    (hidden / "big.bin").write_bytes(b"z" * 50)
    hidden.chmod(0)
    try:
        with pytest.raises(VerifyError) as ei:
            submission_tree_bytes(tmp_path)
        assert ei.value.code is ErrorCode.RESOURCE_EXCEEDED and "notes" in ei.value.message
    finally:
        hidden.chmod(0o755)


def test_tree_bytes_skips_a_file_that_vanished_between_listing_and_stat(tmp_path, monkeypatch):
    # A search loop rewriting a scratch file under submission/ while the verifier walks must not turn
    # into an internal error: a file that is gone by the time it is stat'ed weighs nothing.
    (tmp_path / "keep.ec").write_bytes(b"k" * 10)
    (tmp_path / "tmp.ec").write_bytes(b"t" * 10)
    real = os.lstat

    def racy(p, *a, **k):
        if os.path.basename(p) == "tmp.ec":
            raise FileNotFoundError(2, "gone", p)
        return real(p, *a, **k)
    monkeypatch.setattr(os, "lstat", racy)
    assert submission_tree_bytes(tmp_path) == 10


def test_board_dirs_cap_the_nearest_submission_ancestor_else_the_board_itself(tmp_path, monkeypatch):
    # The platform caps submission/. A board under any directory named submission is capped on that
    # directory; a board elsewhere (a $TMPDIR/board next to sieve caches) is capped on itself alone.
    from ec_verify import files
    from ec_verify.board import verify_board_dirs
    from ec_verify.verify import VerifyConfig
    monkeypatch.setattr(files, "MAX_BOARD_BYTES", 40)
    scratch = tmp_path / "scratch"
    (scratch / "board").mkdir(parents=True)
    (scratch / "sieve-cache.bin").write_bytes(b"s" * 50)
    verify_board_dirs(scratch / "board", None, VerifyConfig())          # siblings of a scratch board are not the submission
    sub = tmp_path / "submission"
    (sub / "boards" / "try3").mkdir(parents=True)
    (sub / "notes.txt").write_bytes(b"n" * 50)
    with pytest.raises(VerifyError) as ei:
        verify_board_dirs(sub / "boards" / "try3", None, VerifyConfig())
    assert ei.value.code is ErrorCode.RESOURCE_EXCEEDED and "submission/ totals 50 bytes" in ei.value.message
