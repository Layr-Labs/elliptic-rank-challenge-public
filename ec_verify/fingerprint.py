"""Curve fingerprint for the novelty gate.

|a_p| for every prime p <= PRIME_BOUND (168 primes), with -1 marking a prime of bad
reduction. A candidate matches a known curve when they agree at every prime that is good
for both and at least MIN_COMPARED primes were compared. Isomorphic curves have equal a_p;
isogenous curves have equal a_p (same L-function); quadratic twists have equal |a_p|. Two
non-isogenous, non-twist curves agreeing on 30+ values of |a_p| do not occur in practice
(each |a_p| takes ~2*sqrt(p) values, so the chance is about prod 1/(4 sqrt p) < 1e-30).

Soundness against twists by a smooth d (every p | d is bad on the twist, so it is skipped):
the verifier rejects candidates with fewer than CANDIDATE_MIN_GOOD good primes below
PRIME_BOUND (CURVE_TOO_MANY_BAD_PRIMES). A class member that passes is good at >= 70
primes, at most 168 - known_good of which are bad for the known curve, so
compared >= 70 - (168 - known_good) >= MIN_COMPARED whenever known_good >= 128.
tests/test_novelty.py asserts that bound over the whole frozen set (minimum 139 in the
2026-10-05 freeze). The floor is a backstop, not a bound we expect to meet: the sharpest
known case, icarm-51 twisted by its 83 smallest good primes (the regression test), still
has 74 good primes; a curve found by search has nowhere near 98 bad primes below 1000.

History. Until 2026-10 the bound was 200 (46 primes) with a fixed 30: high-rank curves have
17+ bad primes below 200, so four of the 710 frozen curves could not reach 30 comparisons
and matched nothing, themselves included. The 2026-10-05 patch went to 1000 and added
"compared >= half of the known curve's good primes"; that is an extra condition on matching
and widened the twist gap it meant to close (icarm-51 twisted by its 83 smallest good
primes was accepted). It was removed the same week; this file is the fixed-floor rule.
"""

from __future__ import annotations

from . import gp as gpmod
from .primes import primes_upto

PRIME_BOUND = 1000
MIN_COMPARED = 30          # agreements needed for a match
CANDIDATE_MIN_GOOD = 70    # good primes below PRIME_BOUND a submission must have (see docstring)
TIMEOUT_S = 60.0


PRIMES = primes_upto(PRIME_BOUND)


def fingerprint(ainvs, gp_bin=None) -> list:
    vec = gpmod.gp_vec([gpmod.gp_int(a) for a in ainvs])
    body = [f"E = ellinit({vec});", "D = E.disc;"]
    for p in PRIMES:
        body.append(f'print("F{p}=", if(D % {p} == 0, -1, abs(ellap(E, {p}))));')
    out = gpmod.run_gp("\n".join(body), timeout=TIMEOUT_S, gp_bin=gp_bin)
    return [int(out[f"F{p}"]) for p in PRIMES]


def good_primes(fp: list) -> int:
    """Number of primes below PRIME_BOUND at which the fingerprinted curve has good reduction."""
    return sum(1 for v in fp if v >= 0)


def matches(candidate: list, known: list) -> bool:
    """True when `candidate` is the known curve or an isogeny / quadratic twist / re-model of it:
    the two agree at every prime good for both and at least MIN_COMPARED primes were compared."""
    compared = 0
    for x, y in zip(candidate, known):
        if x < 0 or y < 0:
            continue
        if x != y:
            return False
        compared += 1
    return compared >= MIN_COMPARED
