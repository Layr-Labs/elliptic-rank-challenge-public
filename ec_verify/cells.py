"""Cells and the contribution formula (docs/verifier-architecture.md §4).

A curve lives in one cell: its torsion group (one of Mazur's 15, keyed as "1", "n" or "nxm"
in PARI's elltors order) and its certified rank. Against a cell baseline b (the record height,
frozen) the contribution is 100·(b − H)/b when H < b and 0 otherwise; a cell with no baseline
(no certified known curve of that torsion reaches the rank) is a breakthrough worth
100 + 100/(1 + H/H_REF): above 100 for every H and strictly decreasing with no clamp, so a
smaller curve in a breakthrough cell always raises the board total. Rank 0 proves nothing and
is never scored: a rank-0 curve contributes 0, never holds a cell, and rank-0 cells are never
built. The total is a pure function of the board files and ec_verify/records.json.
"""

from __future__ import annotations

H_REF = 50.0             # height scale of the breakthrough term: a curve at H_REF earns 150, at 0 earns 200
MAX_BOARD_FILES = 512     # above the 156 known cells plus breakthroughs; the 1200 s budget binds first
BREAKTHROUGH = 100.0
MIN_SCORED_RANK = 1      # cells below this rank contribute 0 and are never held or built

MAZUR_GROUPS = frozenset({(), (2,), (3,), (4,), (5,), (6,), (7,), (8,), (9,), (10,), (12,),
                          (2, 2), (4, 2), (6, 2), (8, 2)})     # the 15 torsion structures over Q


def torsion_key(structure) -> str:
    t = tuple(int(x) for x in structure)
    if t not in MAZUR_GROUPS:
        raise ValueError(f"not a torsion structure of an elliptic curve over Q: {list(t)}")
    if not t:
        return "1"
    return "x".join(str(x) for x in t)


def cell_key(tkey: str, rank: int) -> str:
    return f"{tkey}/{int(rank)}"


def cell_rank(key: str) -> int:
    return int(key.rsplit("/", 1)[1])


def contribution(H: float, baseline: float | None) -> float:
    if baseline is None:
        return BREAKTHROUGH + 100.0 / (1.0 + H / H_REF)
    if H < baseline:
        return 100.0 * (baseline - H) / baseline
    return 0.0


def curve_contribution(certified_rank: int, H: float, baseline: float | None) -> float:
    """What one curve earns in its cell: 0 below MIN_SCORED_RANK (rank 0 proves nothing), else contribution()."""
    return contribution(H, baseline) if certified_rank >= MIN_SCORED_RANK else 0.0


def board_total(contribs: dict) -> float:
    total = 0.0
    for key in sorted(contribs):
        total += contribs[key]
    return total
