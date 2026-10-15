"""The board budget is a hard deadline: every gp stage's timeout is clipped to the time left, so one
slow curve cannot carry the run past the 1200 s the contract promises (CLAUDE.md, submitting §6)."""

import time

import pytest

from conftest import needs_gp
from util import FAST, base5

from ec_verify import ErrorCode, VerifyError, gp as gpmod
from ec_verify.verify import VerifyConfig, verify_text

pytestmark = needs_gp


def test_run_gp_clips_its_timeout_to_the_deadline():
    # Without a deadline a trivial script runs; with one already passed it is a GPTimeout before
    # gp even launches, and with one a few hundred ms away gp is killed at the deadline, not at
    # the stage's own 60 s cap.
    assert gpmod.run_gp("x=1;\nprint(\"X=\", x);\n", timeout=60)["X"] == "1"
    with gpmod.deadline(time.monotonic() - 1.0):
        with pytest.raises(gpmod.GPTimeout, match="budget"):
            gpmod.run_gp("x=1;\n", timeout=60)
    t0 = time.monotonic()
    with gpmod.deadline(t0 + 0.5):
        with pytest.raises(gpmod.GPTimeout, match="budget"):          # names the budget, not a slow stage
            gpmod.run_gp("while(1, x=1);\n", timeout=60)
    assert time.monotonic() - t0 < 10.0
    assert gpmod.run_gp("x=2;\nprint(\"X=\", x);\n", timeout=60)["X"] == "2"     # scope ended: no deadline


def test_verify_text_honours_the_config_deadline():
    with pytest.raises(VerifyError) as ei:
        verify_text(base5(), VerifyConfig(gp_bin=FAST.gp_bin, novelty_enabled=False, conductor_enabled=False,
                                           deadline=time.monotonic() - 1.0))
    assert ei.value.code is ErrorCode.RESOURCE_EXCEEDED
    # the deadline does not leak out of verify_text
    assert gpmod.run_gp("x=3;\nprint(\"X=\", x);\n", timeout=60)["X"] == "3"
