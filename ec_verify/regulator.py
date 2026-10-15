"""Certified independence via the canonical height pairing.

M = ellheightmatrix(E, P_1..P_r) is the Gram matrix of the points under the
Néron-Tate pairing; the points are Z-independent modulo torsion iff M is
positive definite. PARI computes M to a requested precision but offers no
error bounds, so the certificate is:

  1. compute M at P1 = 57 and P2 = 115 decimal digits; require every entry
     to agree within EPS (tripwire for the height-accuracy assumption);
  2. treat EPS as the per-entry error of the P2 matrix; Weyl's inequality
     gives |λ_i(M) − λ_i(M̃)| <= ||M − M̃||_2 <= r·EPS, so
     λ_min(M̃) − r·EPS > MARGIN  ⇒  M positive definite  ⇒  independent;
  3. regulator lower bound = Π max(λ_i(M̃) − r·EPS, 0);
  4. tripwire independent of the matrix: for every point, |h(2P) − 4·h(P)| <=
     EPS·max(1, h(P)) (the canonical height is quadratic), and the working
     precision reported by gp must be the one requested.

The single soundness assumption (documented in docs/verifier-architecture.md)
is that ellheight is accurate to well within EPS at P2 digits; measured
agreement between P1 and P2 on record curves is ~1e-55.
"""

from __future__ import annotations

from . import gp as gpmod
from .gp import gp_point as _pt
from .result import ErrorCode, VerifyError

P1 = 57
P2 = 115
EPS = 1e-40
MARGIN = 1e-30
TIMEOUT_S = 300.0
HINT_MAX_COEFF = 1000
HINT_TIMEOUT_S = 20.0
HINT_MAX_DIGITS = 300      # skip the exact relation check on huge points (ellmul cost)
HINT_MAX_RANK = 16


def _real(s: str) -> float:
    return float(s.replace(" ", ""))


def certify(ainvs, points, gp_bin=None) -> dict:
    """Return the certificate dict or raise POINTS_DEPENDENT / HEIGHT_PRECISION_DISAGREE."""
    r = len(points)
    if r == 0:
        return {"rank": 0, "regulator": 1.0, "regulator_lower_bound": 1.0, "lambda_min": None, "agreement": 0.0}
    vec = gpmod.gp_vec([gpmod.gp_int(a) for a in ainvs])
    pts = gpmod.gp_vec([_pt(x, y) for x, y in points])
    body = "\n".join([
        f"default(realprecision, {P1});",
        f"E = ellinit({vec});",
        f"V = {pts};",
        "M1 = ellheightmatrix(E, V);",
        f"default(realprecision, {P2});",
        f"E = ellinit({vec});",
        "M2 = ellheightmatrix(E, V);",
        'print("AGREE=", vecmax(abs(M1 - M2)));',
        "ev = qfjacobi(M2)[1];",
        'print("LMIN=", vecmin(ev));',
        'print("LMAX=", vecmax(ev));',
        'print("DET=", matdet(M2));',
        f"lb = 1.0; for(i = 1, {r}, lb *= max(ev[i] - {r} * {EPS!r}, 0));",
        'print("LB=", lb);',
        'print("PREC=", default(realprecision));',
        # quadraticity tripwire: h(2P) = 4 h(P) for every point
        f"q = 0.0; for(i = 1, {r}, hp = ellheight(E, V[i]); q = max(q, abs(ellheight(E, ellmul(E, V[i], 2)) - 4 * hp) / max(1, hp)));",
        'print("QUAD=", q);',
        # dependence hint: eigenvector of the smallest eigenvalue, scaled to small integers
        "W = qfjacobi(M2)[2]; k = 1; for(i = 1, #ev, if(ev[i] < ev[k], k = i)); v = W[, k];",
        "m = vecmax(abs(v)); if(m == 0, m = 1); v = v / m;",
        'print("EVEC=", v);',
    ])
    out = gpmod.run_gp(body, timeout=TIMEOUT_S, precision=P1, gp_bin=gp_bin)
    agree = _real(out["AGREE"])
    lmin = _real(out["LMIN"])
    det = _real(out["DET"])
    lb = _real(out["LB"])
    quad = _real(out["QUAD"])
    if int(out["PREC"]) < P2:
        raise VerifyError(ErrorCode.HEIGHT_PRECISION_DISAGREE,
                          f"gp working precision is {out['PREC']} digits, expected {P2}")
    if agree > EPS:
        raise VerifyError(ErrorCode.HEIGHT_PRECISION_DISAGREE,
                          f"height matrices at {P1} and {P2} digits differ by {agree:.3g} > {EPS:g}")
    if quad > EPS:
        raise VerifyError(ErrorCode.HEIGHT_PRECISION_DISAGREE,
                          f"canonical height is not quadratic to working precision: max |h(2P)-4h(P)|/h = {quad:.3g}")
    if lmin - r * EPS <= MARGIN:
        hint = _relation_hint(ainvs, points, out.get("EVEC", ""), gp_bin)
        msg = (f"points are not provably independent: smallest eigenvalue of the height pairing is "
               f"{lmin:.3g} (need > {MARGIN:g} after the {r}·{EPS:g} error allowance); regulator ≈ {det:.3g}")
        if hint:
            msg += f"; verified relation: {hint}"
        raise VerifyError(ErrorCode.POINTS_DEPENDENT, msg)
    return {"rank": r, "regulator": det, "regulator_lower_bound": lb, "lambda_min": lmin,
            "agreement": max(agree, quad)}


def _relation_hint(ainvs, points, evec: str, gp_bin) -> str:
    """Try to turn the near-null eigenvector into an exact integer relation Σ n_i P_i ∈ E(Q)_tors."""
    if len(points) > HINT_MAX_RANK:
        return ""
    if any(max(len(str(c.numerator)), len(str(c.denominator))) > HINT_MAX_DIGITS for x, y in points for c in (x, y)):
        return ""
    try:
        comps = [float(t.replace(" ", "")) for t in evec.strip("[]~ ").split(",")]
    except ValueError:
        return ""
    if not comps or len(comps) != len(points):
        return ""
    nz = [abs(c) for c in comps if abs(c) > 1e-9]
    if not nz:
        return ""
    # One exact check only (bounded cost): the scaling whose rounding error is smallest.
    best = None
    for scale in (1, 2, 3, 4, 6, 12):
        raw = [c / min(nz) * scale for c in comps]
        ns = [round(v) for v in raw]
        if all(n == 0 for n in ns) or max(abs(n) for n in ns) > HINT_MAX_COEFF:
            continue
        err = max(abs(v - n) for v, n in zip(raw, ns))
        if best is None or err < best[0]:
            best = (err, ns)
    if best is None or best[0] > 1e-6:
        return ""
    ns = best[1]
    if _relation_holds(ainvs, points, ns, gp_bin):
        return " + ".join(f"{n}·P{i + 1}" for i, n in enumerate(ns) if n) + " ∈ E(Q)_tors"
    return ""


def _relation_holds(ainvs, points, ns, gp_bin) -> bool:
    vec = gpmod.gp_vec([gpmod.gp_int(a) for a in ainvs])
    body = [f"E = ellinit({vec});", "S = [0];"]
    for (x, y), n in zip(points, ns):
        if n:
            body.append(f"S = elladd(E, S, ellmul(E, {_pt(x, y)}, {gpmod.gp_int(n)}));")
    body.append('print("TORS=", if(S == [0], 1, ellorder(E, S)));')
    try:
        out = gpmod.run_gp("\n".join(body), timeout=HINT_TIMEOUT_S, gp_bin=gp_bin)
        return int(out["TORS"]) > 0
    except (gpmod.GPError, KeyError, ValueError):
        return False
