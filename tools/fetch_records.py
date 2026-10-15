#!/usr/bin/env python3
"""Maintainer tool: parse Andrej Dujella's rank-record pages into fixtures.

Reads cached HTML (default: tools/cache, or --fetch to download
from https://web.math.pmf.unizg.hr/~duje/tors/) and writes

  tests/corpus/<page>.ec       the first curve on the page with its listed
                               independent points, in submission format
                               (torsion points are not listed; the verifier drops them anyway)
  tools/anchor_curves.json     a-invariants + published rank for every anchor
                               curve (classics ranks 0-7 included inline); its
                               LMFDB rows enter the known set via tools/fetch_known.py

Pages: rk{20,24,25,26,27,28,29,30,31}.html (rank >= N records) and rkeq{10,...}.html
(rank = N exactly). The exact-rank pages list two curves; the first is taken here (parse_page).
parse_curves, used by tools/fetch_known.py for the raw mirror, reads every curve on a page,
including the per-torsion pages z*.html.
"""

from __future__ import annotations

import argparse
import html
import itertools
import json
import pathlib
import re
import sys
import urllib.request
from fractions import Fraction

ROOT = pathlib.Path(__file__).resolve().parent.parent
BASE = "https://web.math.pmf.unizg.hr/~duje/tors/"
RECORD_PAGES = {"rk20": 20, "rk24": 24, "rk25": 25, "rk26": 26, "rk27": 27, "rk28": 28, "rk29": 29, "rk30": 30, "rk31": 31}
EXACT_PAGES = {"rkeq10": 10, "rkeq11": 11, "rkeq14": 14, "rkeq15": 15, "rkeq17": 17,
               "rkeq18": 18, "rkeq19": 19, "rkeq20": 20}
CLASSICS = [  # LMFDB labels; ranks 4-7 confirmed with ellrank during the design pass
    ("11a1", 0, [0, -1, 1, -10, -20]),
    ("37a1", 1, [0, 0, 1, -1, 0]),
    ("389a1", 2, [0, 1, 1, -2, 0]),
    ("5077a1", 3, [0, 0, 1, -7, 6]),
    ("234446a1", 4, [1, -1, 0, -79, 289]),
    ("19047851a1", 5, [0, 0, 1, -79, 342]),
    ("5187563742a1", 6, [1, 1, 0, -2582, 48720]),
    ("382623908456a1", 7, [0, 0, 0, -10012, 346900]),
]


def strip_html(raw: str) -> str:
    t = html.unescape(re.sub(r"<[^>]*>", "", raw).replace("&nbsp;", " "))
    return re.sub(r"\s+", " ", t)


def parse_equation(eq: str) -> list[int]:
    """`y2 [+ xy] [+ y] = x3 [± x2] ± A x ± B` -> [a1, a2, a3, a4, a6]."""
    lhs, rhs = eq.split("=", 1)
    a1 = 1 if re.search(r"\+\s*xy", lhs) else 0
    a3 = 1 if re.search(r"\+\s*y\b", lhs.replace("xy", "")) else 0
    rhs = rhs.replace(" ", "")
    m = re.fullmatch(r"x3([+-]x2)?(?:([+-]\d*)x)?([+-]\d+)?", rhs)
    if not m:
        raise ValueError(f"unrecognised right-hand side: {rhs[:80]}")
    a2 = 0 if not m.group(1) else (1 if m.group(1) == "+x2" else -1)
    a4 = 0
    if m.group(2):
        s = m.group(2)
        a4 = int(s if len(s) > 1 else s + "1")
    a6 = int(m.group(3)) if m.group(3) else 0
    return [a1, a2, a3, a4, a6]


_NUM = r"-?\d(?:[\d\s\\]*\d)?(?:\s*/\s*\d(?:[\d\s\\]*\d)?)?"      # a big integer or fraction, possibly broken by a line break or a stray backslash
_POINT = re.compile(r"P\d+\s*=\s*\[\s*(" + _NUM + r")\s*,\s*(" + _NUM + r")\s*\]?")   # one page forgets a closing bracket


def _clean_num(s: str) -> str:
    return re.sub(r"[\s\\]", "", s)


def parse_points(text: str) -> list[tuple[str, str]]:
    """'P1 = [x, y]' pairs; coordinates may be fractions. A number broken by whitespace or a stray
    backslash (both occur on the pages) is re-joined; parse_curves then checks every point on the
    curve, so a wrong join is dropped there rather than mirrored."""
    return [(_clean_num(m.group(1)), _clean_num(m.group(2))) for m in _POINT.finditer(text)]


def parse_page(raw: str, rank: int):
    t = strip_html(raw)
    i = t.find("y2")
    j = t.find("Independent points", i)
    if i < 0 or j < 0:
        raise ValueError("page has no equation / point list")
    eq = t[i:j]
    k = eq.find("Torsion")
    if k >= 0:
        eq = eq[:k]
    ainvs = parse_equation(eq.strip())
    pts = parse_points(t[j:])[:rank]   # first curve only on two-curve pages
    return ainvs, pts


def on_curve(ainvs, x, y) -> bool:
    """Exact check of y² + a1xy + a3y = x³ + a2x² + a4x + a6 for rational x, y ('num/den' allowed)."""
    a1, a2, a3, a4, a6 = (Fraction(int(a)) for a in ainvs)
    X, Y = Fraction(str(x)), Fraction(str(y))
    return Y * Y + a1 * X * Y + a3 * Y == X ** 3 + a2 * X * X + a4 * X + a6


_TERM = re.compile(r"([+-]+)(\d*)(x2|x)?")
_SLOT = {"x2": 0, "x": 1, None: 2}


def rhs_candidates(rhs: str) -> list:
    """'- x2 - -5x + 3' -> every reading of the right-hand side as [a2, a4, a6]: one list normally,
    two per doubled sign (a typo on a few of Dujella's pages; the listed points decide)."""
    s = rhs.replace(" ", "")
    pos, terms = 0, []
    for m in _TERM.finditer(s):
        if m.start() != pos:
            raise ValueError(f"unrecognised right-hand side: {rhs[:80]}")
        pos = m.end()
        signs, digits, var = m.groups()
        if not digits and var is None:
            raise ValueError(f"unrecognised right-hand side: {rhs[:80]}")
        terms.append((signs, int(digits) if digits else 1, var))
    if pos != len(s) or not terms:
        raise ValueError(f"unrecognised right-hand side: {rhs[:80]}")
    if len({v for _, _, v in terms}) != len(terms):
        raise ValueError(f"repeated term in right-hand side: {rhs[:80]}")
    choices = [((1,) if signs == "+" else (-1,) if signs == "-" else (1, -1)) for signs, _, _ in terms]
    out = []
    for picks in itertools.product(*choices):
        a = [0, 0, 0]
        for sign, (_, coef, var) in zip(picks, terms):
            a[_SLOT[var]] = sign * coef
        out.append(a)
    return out


_EQ = re.compile(r"y\s*2\s*((?:\+\s*xy)?)\s*((?:\+\s*y)?)\s*=\s*x\s*3\s*(.*?)\s*(?=Torsion points|Independent points)", re.S)
_MARK = re.compile(r"P\d+\s*=")
_ANY_EQ = re.compile(r"y\s*2[^=]{0,12}=\s*x\s*3")     # every Weierstrass equation, with or without a point list after it


def parse_curves(raw: str, problems: list | None = None) -> list:
    """Every curve on a Dujella page (the per-torsion pages list dozens: each an equation, its torsion
    points and 'Independent points of infinite order: P1 = [..] ...'). Torsion points are not
    returned. An equation without a point list, or a point list that does not parse completely,
    raises: a dropped point understates the published rank. Every listed point is checked on the
    curve exactly. A listed point that is not on the curve (a typo on the page) is dropped and the
    curve kept with `listed` = the announced count; a curve none of whose points fit any reading
    of its equation (the typo is in the equation) is omitted. Both are appended to `problems` as
    dicts (`kind` 'dropped' or 'omitted', `index`, `ainvs` of the best reading, `message`) when a
    list is given and raise otherwise; the caller decides whether an omission is acceptable: a
    mis-parsed coefficient would mirror a curve nobody published."""
    t = strip_html(raw)
    eqs = list(_EQ.finditer(t))
    announced = len(_ANY_EQ.findall(t))
    if len(eqs) != announced:
        raise ValueError(f"{announced} equations on the page, {len(eqs)} followed by a point list")
    out = []
    for i, m in enumerate(eqs):
        end = eqs[i + 1].start() if i + 1 < len(eqs) else len(t)
        block = t[m.end():end]
        j = block.find("Independent points")
        if j < 0:
            raise ValueError(f"curve {i + 1}: no independent point list")
        seg = block[j:]
        pts = parse_points(seg)
        marks = len(_MARK.findall(seg))
        if not pts or len(pts) != marks:
            raise ValueError(f"curve {i + 1}: {marks} points announced, {len(pts)} parsed")
        a1 = 1 if m.group(1) else 0
        a3 = 1 if m.group(2) else 0
        cands = [[a1, a2, a3, a4, a6] for a2, a4, a6 in rhs_candidates(m.group(3))]
        fits = [(c, [p for p in pts if on_curve(c, *p)]) for c in cands]
        fits.sort(key=lambda cf: -len(cf[1]))
        ainvs, good = fits[0]
        head = m.group(0)[:60]
        if not good or (len(fits) > 1 and len(fits[1][1]) == len(good)):
            msg = f"curve {i + 1} ({head}...): listed points are not on the curve for any reading of the equation; omitted"
            if problems is None:
                raise ValueError(msg)
            problems.append({"kind": "omitted", "index": i + 1, "ainvs": ainvs, "message": msg})
            continue
        if len(good) != len(pts):
            msg = f"curve {i + 1} ({head}...): {len(pts) - len(good)} of {len(pts)} listed points are not on the curve; dropped"
            if problems is None:
                raise ValueError(msg)
            problems.append({"kind": "dropped", "index": i + 1, "ainvs": ainvs, "message": msg})
        out.append({"ainvs": ainvs, "points": good, "listed": len(pts)})
    return out


def submission_text(ainvs, pts, comment: str) -> str:
    lines = [f"# {comment}", "E " + " ".join(str(a) for a in ainvs)]
    lines += [f"P {x} {y}" for x, y in pts]
    lines.append(f"#CLAIM rank={len(pts)}")
    return "\n".join(lines) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", default=str(ROOT / "tools" / "cache"))
    ap.add_argument("--fetch", action="store_true", help="download pages instead of reading the cache")
    args = ap.parse_args()
    cache = pathlib.Path(args.cache)
    anchors = [{"id": lab, "rank": r, "ainvs": a, "source": "lmfdb"} for lab, r, a in CLASSICS]
    corpus = ROOT / "tests" / "corpus"
    corpus.mkdir(parents=True, exist_ok=True)
    for page, rank in {**EXACT_PAGES, **RECORD_PAGES}.items():
        if args.fetch:
            raw = urllib.request.urlopen(BASE + page + ".html", timeout=60).read().decode("latin-1")
        else:
            raw = (cache / f"{page}.html").read_text(errors="replace")
        ainvs, pts = parse_page(raw, rank)
        kind = "rank >= " if page in RECORD_PAGES else "rank = "
        anchors.append({"id": page, "rank": rank, "ainvs": ainvs,
                        "source": f"dujella {page}.html ({kind}{rank})"})
        if len(pts) != rank:
            print(f"!! {page}: {len(pts)} points, expected {rank}", file=sys.stderr)
        (corpus / f"{page}.ec").write_text(
            submission_text(ainvs, pts, f"Dujella {page}.html: {kind}{rank} curve with its published points"),
            encoding="ascii")
        print(f"{page}: rank {rank}, a4 {len(str(abs(ainvs[3])))} digits, {len(pts)} points")
    (ROOT / "tools" / "anchor_curves.json").write_text(json.dumps(anchors, indent=1) + "\n")
    print(f"wrote {len(anchors)} anchors -> tools/anchor_curves.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
