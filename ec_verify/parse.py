"""Submission grammar (docs/submitting.md §2).

    E <a1> <a2> <a3> <a4> <a6>     integers (each <= MAX_COEFF_DIGITS digits)
    P <x> <y>                      rationals num[/den] (each component <= MAX_POINT_DIGITS digits)
    #CLAIM rank=<k>                optional, recorded only
    # comment / blank              ignored

Fields are separated by spaces/tabs only; integer tokens are ASCII `-?[0-9]+`.
Anything else is PARSE_SYNTAX.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from fractions import Fraction

from .result import ErrorCode, VerifyError

MAX_COEFF_DIGITS = 400
MAX_POINT_DIGITS = 2000
MAX_POINTS = 64
MAX_LINE_CHARS = 10_000
MAX_LINES = 2_000

_INT = re.compile(r"^-?[0-9]+$")
_RAT = re.compile(r"^(-?[0-9]+)(?:/([0-9]+))?$")
_CLAIM = re.compile(r"^#CLAIM[ \t]+rank[ \t]*=[ \t]*([0-9]{1,6})[ \t]*$")
_SEP = re.compile(r"[ \t]+")


@dataclass(frozen=True)
class Submission:
    ainvs: tuple            # (a1, a2, a3, a4, a6) as ints
    points: tuple           # ((Fraction x, Fraction y), ...) in file order
    claimed_rank: int | None
    point_lines: tuple      # original 1-based line numbers of the P lines


def _int(tok: str, what: str, max_digits: int, code: ErrorCode, line: int) -> int:
    if not _INT.match(tok):
        raise VerifyError(ErrorCode.PARSE_SYNTAX, f"line {line}: {what} is not an integer literal: {tok[:40]!r}")
    if len(tok.lstrip("-")) > max_digits:
        raise VerifyError(code, f"line {line}: {what} exceeds {max_digits} digits")
    return int(tok)


def _rat(tok: str, what: str, line: int) -> Fraction:
    m = _RAT.match(tok)
    if not m:
        raise VerifyError(ErrorCode.PARSE_SYNTAX, f"line {line}: {what} is not a rational literal: {tok[:40]!r}")
    num, den = m.group(1), m.group(2) or "1"
    if len(num.lstrip("-")) > MAX_POINT_DIGITS or len(den) > MAX_POINT_DIGITS:
        raise VerifyError(ErrorCode.POINT_TOO_LARGE, f"line {line}: {what} exceeds {MAX_POINT_DIGITS} digits")
    if den == "0" or (len(den) > 1 and den.startswith("0")):
        raise VerifyError(ErrorCode.PARSE_SYNTAX, f"line {line}: {what} has an invalid denominator")
    return Fraction(int(num), int(den))


def parse(text: str) -> Submission:
    if not text.isascii():
        raise VerifyError(ErrorCode.PARSE_SYNTAX, "submission is not ASCII")
    lines = text.split("\n")
    if len(lines) > MAX_LINES:
        raise VerifyError(ErrorCode.PARSE_SYNTAX, f"more than {MAX_LINES} lines")
    ainvs = None
    points: list = []
    point_lines: list = []
    claimed = None
    for no, raw in enumerate(lines, 1):
        if len(raw) > MAX_LINE_CHARS:
            raise VerifyError(ErrorCode.PARSE_SYNTAX, f"line {no}: longer than {MAX_LINE_CHARS} characters")
        line = raw.rstrip("\r")
        if not line.strip(" \t"):
            continue
        if line.startswith("#"):
            if line.startswith("#CLAIM"):
                m = _CLAIM.match(line)
                if not m:
                    raise VerifyError(ErrorCode.PARSE_SYNTAX, f"line {no}: malformed #CLAIM (expected `#CLAIM rank=<1-6 digits>`)")
                if claimed is not None:
                    raise VerifyError(ErrorCode.PARSE_SYNTAX, f"line {no}: duplicate #CLAIM")
                claimed = int(m.group(1))
            continue
        toks = [t for t in _SEP.split(line.strip(" \t")) if t]
        if toks[0] == "E":
            if ainvs is not None:
                raise VerifyError(ErrorCode.PARSE_SYNTAX, f"line {no}: second E line")
            if len(toks) != 6:
                raise VerifyError(ErrorCode.PARSE_SYNTAX, f"line {no}: E needs exactly five a-invariants")
            ainvs = tuple(_int(t, f"a-invariant {i}", MAX_COEFF_DIGITS, ErrorCode.CURVE_TOO_LARGE, no)
                          for i, t in enumerate(toks[1:], 1))
        elif toks[0] == "P":
            if ainvs is None:
                raise VerifyError(ErrorCode.PARSE_SYNTAX, f"line {no}: P before E")
            if len(toks) != 3:
                raise VerifyError(ErrorCode.PARSE_SYNTAX, f"line {no}: P needs exactly x and y")
            if len(points) >= MAX_POINTS:
                raise VerifyError(ErrorCode.TOO_MANY_POINTS, f"more than {MAX_POINTS} points")
            points.append((_rat(toks[1], "x", no), _rat(toks[2], "y", no)))
            point_lines.append(no)
        else:
            raise VerifyError(ErrorCode.PARSE_SYNTAX, f"line {no}: unknown record {toks[0][:20]!r}")
    if ainvs is None:
        raise VerifyError(ErrorCode.PARSE_EMPTY, "no E line")
    return Submission(ainvs=ainvs, points=tuple(points), claimed_rank=claimed, point_lines=tuple(point_lines))
