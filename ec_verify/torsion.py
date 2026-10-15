"""Torsion: exact order of E(Q)_tors and exact stripping of finite-order points."""

from __future__ import annotations

from . import gp as gpmod
from .gp import gp_point as _pt

TIMEOUT_S = 120.0


def parse_gp_ints(s: str) -> list[int]:
    """'[6, 2]' -> [6, 2]; '[]' -> []."""
    body = s.strip()
    if not (body.startswith("[") and body.endswith("]")):
        raise gpmod.GPError(f"GP_ERROR: expected a gp vector, got {s!r}")
    body = body[1:-1].strip()
    if not body:
        return []
    try:
        return [int(tok) for tok in body.split(",")]      # int("") raises: an empty entry fails closed
    except ValueError as e:
        raise gpmod.GPError(f"GP_ERROR: malformed gp vector {s!r}: {e}")


def strip_torsion(ainvs, points, gp_bin=None):
    """Return (torsion_order, structure, [points of infinite order], [indices dropped]).

    structure is elltors(E)[2]: [] (trivial), [n] (Z/n) or [n, m] (Z/n x Z/m, m | n)."""
    vec = gpmod.gp_vec([gpmod.gp_int(a) for a in ainvs])
    body = [f"E = ellinit({vec});", "T = elltors(E);", 'print("TORS=", T[1]);', 'print("STRUCT=", T[2]);']
    for i, (x, y) in enumerate(points):
        body.append(f'print("O{i}=", ellorder(E, {_pt(x, y)}));')
    out = gpmod.run_gp("\n".join(body), timeout=TIMEOUT_S, gp_bin=gp_bin)
    tors = int(out["TORS"])
    structure = parse_gp_ints(out["STRUCT"])
    keep, dropped = [], []
    for i, p in enumerate(points):
        if int(out[f"O{i}"]) == 0:
            keep.append(p)
        else:
            dropped.append(i)
    return tors, structure, keep, dropped
