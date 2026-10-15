"""Exact point checks (pure Python) and the exact curve invariants."""

from fractions import Fraction

import pytest

from util import corpus, gens_5077a

from ec_verify import ErrorCode, VerifyError
from ec_verify.curve import discriminant, naive_height_log10, require_nonsingular
from ec_verify.parse import parse
from ec_verify.points import check_points, negate, on_curve


def run(text):
    s = parse(text)
    check_points(s.ainvs, s.points, s.point_lines)
    return s


def expect(text, code):
    with pytest.raises(VerifyError) as ei:
        run(text)
    assert ei.value.code is code, ei.value


def test_generators_pass():
    run(gens_5077a())


def test_off_curve():
    expect(gens_5077a().replace("P -2 3", "P -2 4"), ErrorCode.POINT_OFF_CURVE)


def test_duplicate_point():
    expect(gens_5077a().replace("P -1 3", "P -2 3"), ErrorCode.POINT_DUPLICATE)


def test_negative_of_point():
    # -(x, y) = (x, -y - a1 x - a3): on 5077a, -(-2, 3) = (-2, -4)
    assert negate((0, 0, 1, -7, 6), Fraction(-2), Fraction(3)) == (Fraction(-2), Fraction(-4))
    assert on_curve((0, 0, 1, -7, 6), Fraction(-2), Fraction(-4))
    expect(gens_5077a().replace("P -1 3", "P -2 -4"), ErrorCode.POINT_DUPLICATE)


def test_rational_points_on_record_curve():
    s = run(corpus("rk30.ec"))
    assert any(p[0].denominator > 1 for p in s.points), "rk30 has points with denominators"
    assert len(s.points) == 30


def test_discriminant_known_values():
    assert discriminant((0, 0, 1, -7, 6)) == 5077          # 5077a1
    assert discriminant((0, -1, 1, -10, -20)) == -161051   # 11a1 = -11^5
    assert discriminant((0, 0, 1, -1, 0)) == 37            # 37a1


def test_singular_rejected():
    with pytest.raises(VerifyError) as ei:
        require_nonsingular((0, 0, 0, 0, 0))
    assert ei.value.code is ErrorCode.CURVE_SINGULAR
    with pytest.raises(VerifyError):
        require_nonsingular((0, 0, 0, -3, 2))   # y^2 = (x-1)^2 (x+2)


def test_naive_height():
    assert naive_height_log10((0, 0, 1, -7, 6)) == pytest.approx(0.845098, abs=1e-6)
