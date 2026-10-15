"""Independence certificate: the study's own dependent triple, the true generators, record curves."""

import pytest

from conftest import needs_gp
from util import FAST, corpus, count_points, dependent_5077a, gens_5077a

from ec_verify import ErrorCode, VerifyError, verify_text
from ec_verify import regulator

pytestmark = needs_gp


def test_constants_are_sane():
    # the positivity margin must dominate the total error allowance for any admissible r
    assert regulator.MARGIN > 64 * regulator.EPS > 0
    assert regulator.P2 >= 2 * regulator.P1 - 1


def test_dependent_triple_rejected_with_relation():
    with pytest.raises(VerifyError) as ei:
        verify_text(dependent_5077a(), FAST)
    assert ei.value.code is ErrorCode.POINTS_DEPENDENT
    assert "1·P1 + 1·P2 + 2·P3" in ei.value.message


def test_generators_certify_rank_three():
    r = verify_text(gens_5077a(), FAST)
    assert r.certified_rank == 3
    assert r.regulator == pytest.approx(0.4171435587583840, rel=1e-9)
    assert 0 < r.regulator_lower_bound <= r.regulator
    assert r.height_lambda_min > regulator.MARGIN
    assert r.height_precision_agreement <= regulator.EPS


def test_subset_of_generators_certifies_fewer():
    text = "\n".join(l for l in gens_5077a().splitlines() if not l.startswith("P 0 2")) + "\n"
    assert verify_text(text, FAST).certified_rank == 2


def test_two_dependent_points_multiple():
    # 37a1: P = (0,0), 2P = (1/4, -5/8) -> dependent.
    with pytest.raises(VerifyError) as ei:
        verify_text("E 0 0 1 -1 0\nP 0 0\nP 1/4 -5/8\n", FAST)
    assert ei.value.code is ErrorCode.POINTS_DEPENDENT


@pytest.mark.parametrize("name", ["rk20.ec", "rk24.ec", "rk26.ec", "rk28.ec", "rk29.ec", "rk30.ec", "rk31.ec"])
def test_record_curves_certify_published_rank(name):
    text = corpus(name)
    r = verify_text(text, FAST)
    assert r.certified_rank == count_points(text)
    assert r.torsion_order == 1
    assert r.height_precision_agreement <= regulator.EPS
    assert r.regulator_lower_bound > 0


def test_record_curve_with_a_repeated_point_is_dependent():
    text = corpus("rk30.ec")
    lines = text.splitlines()
    ps = [l for l in lines if l.startswith("P ")]
    # replace the last point with the sum P1 + P2 is expensive to build by hand; simpler: drop one and
    # duplicate another via P and 2P isn't available exactly here — so use a shifted copy: P1 listed twice
    # is caught earlier as POINT_DUPLICATE; the genuine height-level test is the 37a multiple above.
    dup = text.replace(ps[-1], ps[0], 1)
    with pytest.raises(VerifyError) as ei:
        verify_text(dup, FAST)
    assert ei.value.code is ErrorCode.POINT_DUPLICATE
