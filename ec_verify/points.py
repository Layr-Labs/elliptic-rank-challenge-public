"""Exact rational-point checks (pure Python, no PARI)."""

from __future__ import annotations

from fractions import Fraction

from .result import ErrorCode, VerifyError


def on_curve(ainvs, x: Fraction, y: Fraction) -> bool:
    a1, a2, a3, a4, a6 = ainvs
    return y * y + a1 * x * y + a3 * y == x ** 3 + a2 * x * x + a4 * x + a6


def negate(ainvs, x: Fraction, y: Fraction):
    a1, _, a3, _, _ = ainvs
    return (x, -y - a1 * x - a3)


def check_points(ainvs, points, point_lines):
    """Every point on the curve; no P repeated or listed with -P."""
    seen: dict = {}
    for (x, y), no in zip(points, point_lines):
        if not on_curve(ainvs, x, y):
            raise VerifyError(ErrorCode.POINT_OFF_CURVE, f"line {no}: point does not satisfy the curve equation")
        key = (x, y)
        nkey = negate(ainvs, x, y)
        if key in seen:
            raise VerifyError(ErrorCode.POINT_DUPLICATE, f"line {no}: repeats the point on line {seen[key]}")
        if nkey in seen:
            raise VerifyError(ErrorCode.POINT_DUPLICATE, f"line {no}: is the negative of the point on line {seen[nkey]}")
        seen[key] = no
