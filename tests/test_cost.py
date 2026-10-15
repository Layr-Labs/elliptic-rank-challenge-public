"""Cost bounds: a hostile-but-legal file must be rejected quickly, not run the clock."""

import time

import pytest

from conftest import needs_gp
from util import FAST

from ec_verify import ErrorCode, VerifyError, verify_text
from ec_verify import gp as gpmod

pytestmark = [needs_gp, pytest.mark.slow]


def test_64_dependent_huge_points_reject_within_budget():
    body = ["E = ellinit([0,0,1,-1,0]);", "P = [0,0];"]
    ks = list(range(150, 214))
    for k in ks:
        body.append(f'Q = ellmul(E, P, {k}); print("X{k}=", numerator(Q[1]), "/", denominator(Q[1])); print("Y{k}=", numerator(Q[2]), "/", denominator(Q[2]));')
    out = gpmod.run_gp("\n".join(body), timeout=300)
    text = "E 0 0 1 -1 0\n" + "".join(f"P {out[f'X{k}'].replace(' ', '')} {out[f'Y{k}'].replace(' ', '')}\n" for k in ks)
    t0 = time.time()
    with pytest.raises(VerifyError) as ei:
        verify_text(text, FAST)
    assert ei.value.code is ErrorCode.POINTS_DEPENDENT
    assert time.time() - t0 < 60
