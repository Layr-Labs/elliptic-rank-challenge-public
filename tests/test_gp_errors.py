"""gp 2.17.4 keeps executing a script after a recoverable error (exit 0, diagnostics on stderr), so the
completion sentinel alone proves nothing. run_gp must turn a gp error into GPError, and a warning
(stack growth, user warning) must not."""

import re
import subprocess

import pytest

from conftest import needs_gp

from ec_verify import gp as gpmod

pytestmark = needs_gp


def test_a_gp_error_mid_script_raises_gperror_not_a_partial_dict():
    with pytest.raises(gpmod.GPError, match="impossible inverse"):
        gpmod.run_gp('x = 1/0;\nprint("AFTER=1");\ny = x + 1;\nprint("Y=", y);\n', timeout=30)


def test_a_gp_warning_is_not_an_error():
    out = gpmod.run_gp('warning("just a warning");\nprint("X=", 1);\n', timeout=30)
    assert out["X"] == "1"


def test_an_error_whose_message_contains_the_word_warning_is_still_an_error():
    # Only gp's own "Warning:" / "user warning:" lines are benign; an error ("user error: ...",
    # or any "***" line) that happens to mention the word must still fail the stage.
    with pytest.raises(gpmod.GPError, match="user error"):
        gpmod.run_gp('error("a warning-shaped message");\nprint("X=", 1);\n', timeout=30)


def _gp_stderr(script: str, tmp_path, parisize: int) -> subprocess.CompletedProcess:
    p = tmp_path / "probe.gp"
    p.write_text(script + 'print("EC_DONE=1");\n')
    cmd = gpmod.gp_command(gpmod.find_gp(), p)
    cmd[cmd.index(f"parisize={gpmod.PARISIZE}")] = f"parisize={parisize}"
    return subprocess.run(cmd, capture_output=True, text=True, timeout=60)


GROWTH = 'v = vector(300000, i, i^20);\nprint("X=", 1);\n'


def test_a_stack_growth_warning_inside_a_pari_function_is_not_an_error(tmp_path, monkeypatch):
    # PARI prefixes a warning raised inside a function with the function's name:
    # "  *** _^s: Warning: increasing stack size to 6400000." A large curve whose height or
    # minimality stage grows the stack inside ellheight/ellminimalmodel must not be GP_ERROR.
    probe = _gp_stderr(GROWTH, tmp_path, parisize=400000)            # the shape really is emitted at this size
    assert probe.returncode == 0 and re.search(r"^\s*\*\*\* \S+: Warning: increasing stack size", probe.stderr, re.M), probe.stderr
    monkeypatch.setattr(gpmod, "PARISIZE", 400000)
    assert gpmod.run_gp(GROWTH, timeout=60)["X"] == "1"


def test_an_error_inside_a_pari_function_is_still_an_error():
    # the same "*** <function>:" prefix on an error ("  *** _/_: impossible inverse") stays fatal,
    # classified by the stderr pass (not merely caught by the missing sentinel)
    with pytest.raises(gpmod.GPError, match=r"^GP_ERROR: gp reported an error: .*impossible inverse"):
        gpmod.run_gp('f(x) = 1/x;\ny = f(0);\nprint("X=", 1);\n', timeout=30)


def test_recover_off_makes_an_error_abort_gp_without_the_sentinel(tmp_path):
    # Structural backstop: even if the stderr classifier missed an error line, gp exits 1 and never
    # prints the sentinel, so run_gp cannot return a partial result.
    probe = _gp_stderr('x = 1/0;\nprint("X=", 1);\n', tmp_path, parisize=gpmod.PARISIZE)
    assert probe.returncode != 0 and "EC_DONE=1" not in probe.stdout
    probe = _gp_stderr('warning("w");\nprint("X=", 1);\n', tmp_path, parisize=gpmod.PARISIZE)
    assert probe.returncode == 0 and "EC_DONE=1" in probe.stdout
