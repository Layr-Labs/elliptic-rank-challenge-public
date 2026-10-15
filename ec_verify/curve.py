"""Exact curve invariants (pure Python) and the small-prime minimality gate."""

from __future__ import annotations

import math

from . import gp as gpmod
from .primes import primes_upto
from .result import ErrorCode, VerifyError

SMALL_PRIME_BOUND = 10 ** 6
LOCALRED_TIMEOUT_S = 120.0
CONDUCTOR_TIMEOUT_S = 20.0
CONDUCTOR_MAX_DISC_DIGITS = 40   # deterministic gate: factoring 40 digits is always fast


def invariants(ainvs):
    a1, a2, a3, a4, a6 = ainvs
    b2 = a1 * a1 + 4 * a2
    b4 = 2 * a4 + a1 * a3
    b6 = a3 * a3 + 4 * a6
    b8 = a1 * a1 * a6 + 4 * a2 * a6 - a1 * a3 * a4 + a2 * a3 * a3 - a4 * a4
    disc = -b2 * b2 * b8 - 8 * b4 ** 3 - 27 * b6 * b6 + 9 * b2 * b4 * b6
    return b2, b4, b6, b8, disc


def discriminant(ainvs) -> int:
    return invariants(ainvs)[4]


def c_invariants(ainvs):
    b2, b4, b6, b8, _ = invariants(ainvs)
    c4 = b2 * b2 - 24 * b4
    c6 = -b2 ** 3 + 36 * b2 * b4 - 216 * b6
    return c4, c6


def naive_height(ainvs) -> float:
    """H = (1/12) log max(|c4|^3, c6^2): invariant under isomorphism (u = ±1 on a model that is
    minimal, as the small-prime gate enforces), so a re-model cannot shrink it."""
    c4, c6 = c_invariants(ainvs)
    m = max(abs(c4) ** 3, c6 * c6)
    if m == 0:
        return 0.0
    # exact-ish log of a huge integer: log(m) = log(mantissa) + digits*log(10)
    s = str(m)
    k = min(len(s), 20)
    return round((math.log(int(s[:k])) + (len(s) - k) * math.log(10)) / 12, 6)


def require_nonsingular(ainvs) -> int:
    d = discriminant(ainvs)
    if d == 0:
        raise VerifyError(ErrorCode.CURVE_SINGULAR, "discriminant is zero; not an elliptic curve")
    return d


def naive_height_log10(ainvs) -> float:
    m = max(abs(a) for a in ainvs)
    return round(math.log10(m), 6) if m > 0 else 0.0


def digits(n: int) -> int:
    return len(str(abs(n)))


_PRIMES = None


def small_bad_primes(disc: int, bound: int = SMALL_PRIME_BOUND):
    global _PRIMES
    if _PRIMES is None:
        _PRIMES = primes_upto(bound)
    d = abs(disc)
    return [p for p in _PRIMES if d % p == 0]


def require_minimal_at_small_primes(ainvs, disc: int, gp_bin=None):
    """Reject a model that is non-minimal at any prime p <= 10^6 dividing Δ.

    Uses elllocalred per prime (Tate's algorithm), never a factorisation of Δ.
    Non-minimality at a larger prime is undetectable without factoring and is
    harmless: heights are model-invariant and ellap reduces locally.
    """
    bad = small_bad_primes(disc)
    if not bad:
        return []
    vec = gpmod.gp_vec([gpmod.gp_int(a) for a in ainvs])
    body = [f"E = ellinit({vec});"]
    for p in bad:
        body.append(f'L = elllocalred(E, {p}); print("U{p}=", L[3][1]);')
    out = gpmod.run_gp("\n".join(body), timeout=LOCALRED_TIMEOUT_S, gp_bin=gp_bin)
    for p in bad:
        u = out.get(f"U{p}")
        if u is None:
            raise gpmod.GPError(f"GP_ERROR: missing elllocalred at {p}")
        if u != "1":
            raise VerifyError(ErrorCode.CURVE_NON_MINIMAL,
                              f"model is not minimal at p={p} (scale factor u={u}); submit the minimal model")
    return bad


def conductor_best_effort(ainvs, disc: int, gp_bin=None):
    """Informational only. Attempted only when |Δ| has <= CONDUCTOR_MAX_DISC_DIGITS digits
    (a deterministic gate); the time cap is a safety net. Returned as a decimal string."""
    if digits(disc) > CONDUCTOR_MAX_DISC_DIGITS:
        return None
    vec = gpmod.gp_vec([gpmod.gp_int(a) for a in ainvs])
    try:
        out = gpmod.run_gp(f"E = ellinit({vec});\nprint(\"N=\", ellglobalred(E)[1]);",
                           timeout=CONDUCTOR_TIMEOUT_S, gp_bin=gp_bin)
        return str(int(out["N"]))
    except (gpmod.GPError, KeyError, ValueError):
        return None
