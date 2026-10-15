"""benchmark.sh's pre-sandbox failures must speak the documented interface: a `REJECTED: <CODE>: <detail>`
line on stdout, exit 1, no score.json. Run on a scratch copy with a stub venv and a stub gp (the
verify stage itself cannot run nested inside a sandbox; these paths exit before it)."""

from __future__ import annotations

import os
import pathlib
import shutil
import subprocess

import pytest

from tests.util import ROOT

pytestmark = pytest.mark.skipif(shutil.which("git") is None, reason="git not available")

_ENV = {**os.environ, "GIT_CONFIG_GLOBAL": os.devnull, "GIT_CONFIG_NOSYSTEM": "1",
        "GIT_AUTHOR_NAME": "t", "GIT_AUTHOR_EMAIL": "t@t", "GIT_COMMITTER_NAME": "t", "GIT_COMMITTER_EMAIL": "t@t"}


def _scratch_repo(tmp_path: pathlib.Path) -> pathlib.Path:
    repo = tmp_path / "repo"
    (repo / "submission" / "board").mkdir(parents=True)
    shutil.copy(ROOT / "benchmark.sh", repo / "benchmark.sh")
    (repo / "tools").mkdir()
    shutil.copy(ROOT / "tools" / "export_parent_board.sh", repo / "tools" / "export_parent_board.sh")
    vpy = repo / ".venv-bench" / "bin" / "python"
    vpy.parent.mkdir(parents=True)
    vpy.write_text("#!/bin/sh\nexit 0\n", encoding="ascii"); vpy.chmod(0o755)
    gp = repo / "tools" / "pari" / "bin" / "gp"
    gp.parent.mkdir(parents=True)
    gp.write_text("#!/bin/sh\necho 2.17.4\n", encoding="ascii"); gp.chmod(0o755)
    (repo / "score.json").write_text("{\"score\": 99}\n", encoding="ascii")     # stale: must be gone after
    subprocess.run(["git", "-C", str(repo), "init", "-q"], check=True, env=_ENV, capture_output=True)
    subprocess.run(["git", "-C", str(repo), "add", "-A"], check=True, env=_ENV, capture_output=True)
    subprocess.run(["git", "-C", str(repo), "commit", "-q", "-m", "init"], check=True, env=_ENV, capture_output=True)
    return repo


def _run(repo: pathlib.Path, **env):
    e = {k: v for k, v in {**_ENV, **env}.items() if v is not None}
    return subprocess.run(["bash", str(repo / "benchmark.sh")], cwd=repo, env=e, capture_output=True, text=True, timeout=60)


def test_bad_parent_ref_is_a_rejected_line(tmp_path):
    repo = _scratch_repo(tmp_path)
    proc = _run(repo, EC_PARENT_REF="HEAD^")                      # the root commit has no parent
    assert proc.returncode == 1
    assert proc.stdout.startswith("REJECTED: BOARD_REGRESSION:") and "HEAD^" in proc.stdout, proc.stdout + proc.stderr
    assert not (repo / "score.json").exists()
