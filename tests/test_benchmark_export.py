"""tools/export_parent_board.sh: the parent-board export is attribute-free (no git archive).

The export runs in benchmark.sh outside the sandbox; it is tested directly on a scratch repo
because benchmark.sh itself cannot run nested inside a sandbox. No gp needed."""

from __future__ import annotations

import os
import pathlib
import shutil
import subprocess

import pytest

from tests.util import ROOT, base5

SCRIPT = ROOT / "tools" / "export_parent_board.sh"

pytestmark = pytest.mark.skipif(shutil.which("git") is None, reason="git not available")

_ENV = {
    **os.environ,
    "GIT_CONFIG_GLOBAL": os.devnull,
    "GIT_CONFIG_NOSYSTEM": "1",
    "GIT_AUTHOR_NAME": "t", "GIT_AUTHOR_EMAIL": "t@t",
    "GIT_COMMITTER_NAME": "t", "GIT_COMMITTER_EMAIL": "t@t",
}


def _git(repo: pathlib.Path, *args: str) -> None:
    subprocess.run(["git", "-C", str(repo), *args], check=True, env=_ENV, capture_output=True)


def _repo(tmp_path: pathlib.Path, with_board: bool = True, symlink: bool = False) -> pathlib.Path:
    repo = tmp_path / "repo"
    sub = repo / "submission"
    board = sub / "board"
    board.mkdir(parents=True)
    _git(repo, "init", "-q")
    _git(repo, "config", "commit.gpgsign", "false")
    (sub / ".gitattributes").write_text("board/*.ec export-ignore\nboard/* export-subst\n", encoding="ascii")
    if with_board:
        (board / "a.ec").write_text(base5(), encoding="ascii")
        if symlink:
            (board / "b.ec").symlink_to("a.ec")
        nested = board / "nested"
        nested.mkdir()
        (nested / "c.ec").write_text(base5(), encoding="ascii")
    else:
        (sub / "README").write_text("no board\n", encoding="ascii")
    _git(repo, "add", "-A")
    _git(repo, "commit", "-q", "-m", "c")
    return repo


def _export(repo: pathlib.Path, dest: pathlib.Path, ref: str = "HEAD") -> subprocess.CompletedProcess:
    return subprocess.run(["bash", str(SCRIPT), str(repo), ref, str(dest)],
                          env=_ENV, capture_output=True, text=True)


def test_export_ignores_gitattributes(tmp_path):
    repo = _repo(tmp_path)
    dest = tmp_path / "out"
    r = _export(repo, dest)
    assert r.returncode == 0, r.stderr
    assert sorted(p.name for p in dest.iterdir()) == ["a.ec"]
    assert (dest / "a.ec").read_bytes() == base5().encode("ascii")
    assert "!!" in r.stderr and "nested/c.ec" in r.stderr


def test_export_skips_symlink_with_note(tmp_path):
    repo = _repo(tmp_path, symlink=True)
    dest = tmp_path / "out"
    r = _export(repo, dest)
    assert r.returncode == 0, r.stderr
    assert sorted(p.name for p in dest.iterdir()) == ["a.ec"]
    assert not (dest / "a.ec").is_symlink()
    assert "!!" in r.stderr and "b.ec" in r.stderr and "120000" in r.stderr


def test_export_without_board_is_empty_dir(tmp_path):
    repo = _repo(tmp_path, with_board=False)
    dest = tmp_path / "out"
    r = _export(repo, dest)
    assert r.returncode == 0, r.stderr
    assert dest.is_dir() and list(dest.iterdir()) == []


def test_export_bad_ref_fails(tmp_path):
    repo = _repo(tmp_path)
    r = _export(repo, tmp_path / "out", ref="deadbeef")
    assert r.returncode == 1
    assert "!!" in r.stderr


def test_benchmark_sh_uses_the_script():
    src = (ROOT / "benchmark.sh").read_text()
    assert "tools/export_parent_board.sh" in src
    assert "git archive" not in src and "archive --format" not in src


def test_export_counts_files_not_lines(tmp_path):
    # A file name with an embedded newline is a legal blob name; it must not break the count
    # (a line count would see two entries and strand every later submission).
    repo = tmp_path / "repo"
    board = repo / "submission" / "board"
    board.mkdir(parents=True)
    _git(repo, "init", "-q")
    (board / "a.ec").write_text(base5(), encoding="ascii")
    pathlib.Path(board, "foo\nbar").write_text("note\n", encoding="ascii")
    _git(repo, "add", "-A")
    _git(repo, "-c", "commit.gpgsign=false", "commit", "-q", "-m", "c")
    dest = tmp_path / "out"
    r = _export(repo, dest)
    assert r.returncode == 0, r.stderr
    assert sorted(p.name for p in dest.iterdir()) == ["a.ec", "foo\nbar"]
    assert (dest / "foo\nbar").read_text() == "note\n"
    assert [p.name for p in tmp_path.iterdir() if "ls-tree" in p.name] == []   # listing removed
