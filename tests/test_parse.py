"""Grammar: every malformed input maps to a stable code; valid input round-trips."""

from fractions import Fraction

import pytest

from util import gens_5077a

from ec_verify import ErrorCode, VerifyError
from ec_verify.parse import MAX_COEFF_DIGITS, MAX_POINTS, parse


def expect(text, code):
    with pytest.raises(VerifyError) as ei:
        parse(text)
    assert ei.value.code is code, f"expected {code}, got {ei.value.code}: {ei.value.message}"


def test_valid_roundtrip():
    s = parse(gens_5077a())
    assert s.ainvs == (0, 0, 1, -7, 6)
    assert s.points == ((Fraction(-2), Fraction(3)), (Fraction(-1), Fraction(3)), (Fraction(0), Fraction(2)))
    assert s.claimed_rank == 3
    assert s.point_lines == (3, 4, 5)


def test_rational_points_and_omitted_denominator():
    s = parse("E 1 0 0 -1 0\nP 1/2 -3/4\nP 7 -1\n")
    assert s.points[0] == (Fraction(1, 2), Fraction(-3, 4))
    assert s.points[1] == (Fraction(7), Fraction(-1))


def test_crlf_and_blank_lines_tolerated():
    s = parse("E 0 0 1 -7 6\r\n\r\nP -2 3\r\n")
    assert len(s.points) == 1


def test_no_claim_is_none():
    assert parse("E 0 0 1 -7 6\n").claimed_rank is None


def test_empty():
    expect("", ErrorCode.PARSE_EMPTY)
    expect("# only a comment\n", ErrorCode.PARSE_EMPTY)


@pytest.mark.parametrize("bad", ["E 0 0 1 -7 +6\n", "E 0 0 1 -7 6.0\n", "E 0 0 1 -7 ６\n", "E 0 0 1 -7 1_0\n",
                                 "E 0 0 1 -7\n", "E 0 0 1 -7 6 0\n", "E 0 0 1 -7 6\nE 0 0 1 -7 6\n",
                                 "P -2 3\nE 0 0 1 -7 6\n", "E 0 0 1 -7 6\nP -2\n", "E 0 0 1 -7 6\nP -2 3 1\n",
                                 "E 0 0 1 -7 6\nP 1/0 3\n", "E 0 0 1 -7 6\nP 1/02 3\n", "E 0 0 1 -7 6\nP 1/-2 3\n",
                                 "E 0 0 1 -7 6\nQ -2 3\n", "E 0 0 1 -7 6\n#CLAIM rank=3\n#CLAIM rank=4\n"])
def test_syntax_errors(bad):
    expect(bad, ErrorCode.PARSE_SYNTAX)


def test_coefficient_cap():
    big = "1" * (MAX_COEFF_DIGITS + 1)
    expect(f"E 0 0 1 -7 {big}\n", ErrorCode.CURVE_TOO_LARGE)
    parse(f"E 0 0 1 -7 {'1' * MAX_COEFF_DIGITS}\n")


def test_point_cap():
    big = "1" * 2001
    expect(f"E 0 0 1 -7 6\nP {big} 3\n", ErrorCode.POINT_TOO_LARGE)
    expect(f"E 0 0 1 -7 6\nP 1/{big} 3\n", ErrorCode.POINT_TOO_LARGE)


def test_too_many_points():
    text = "E 0 0 1 -7 6\n" + "P 0 2\n" * (MAX_POINTS + 1)
    expect(text, ErrorCode.TOO_MANY_POINTS)


def test_claim_too_long_is_syntax_not_crash():
    expect("E 0 0 1 -7 6\n#CLAIM rank=" + "9" * 5000 + "\n", ErrorCode.PARSE_SYNTAX)
    expect("E 0 0 1 -7 6\n#CLAIM rank=1234567\n", ErrorCode.PARSE_SYNTAX)
    expect("E 0 0 1 -7 6\n#CLAIM rank=abc\n", ErrorCode.PARSE_SYNTAX)
    assert parse("E 0 0 1 -7 6\n#CLAIM rank=999999\n").claimed_rank == 999999


@pytest.mark.parametrize("sep", ["\x0b", "\x0c"])
def test_only_space_and_tab_separate_fields(sep):
    expect(f"E{sep}0 0 1 -7 6\n", ErrorCode.PARSE_SYNTAX)
    expect(f"E 0 0 1 -7 6\nP -2{sep}3\n", ErrorCode.PARSE_SYNTAX)
