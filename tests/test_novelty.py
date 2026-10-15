"""The novelty gate: known curves, their isogenies, twists and re-models are rejected; new curves pass."""

import pytest

from conftest import needs_gp
from util import FAST_NOVELTY, base5, corpus, gens_5077a

from ec_verify import ErrorCode, VerifyError, verify_text
from ec_verify import fingerprint as fp_mod, gp as gpmod, known_curves

pytestmark = needs_gp


def expect_known(text):
    with pytest.raises(VerifyError) as ei:
        verify_text(text, FAST_NOVELTY)
    assert ei.value.code is ErrorCode.KNOWN_CURVE, ei.value
    return ei.value.message


def test_frozen_set_is_large_and_covers_the_record():
    t = known_curves.load()
    assert len(t["curves"]) >= 4000
    assert max(c["rank_claimed"] for c in t["curves"]) == 31
    assert t["primes_upto"] == fp_mod.PRIME_BOUND == 1000
    assert all(len(c["fingerprint"]) == len(fp_mod.PRIMES) for c in t["curves"])


def test_every_frozen_curve_matches_its_own_fingerprint():
    """The gate must catch every curve it is supposed to know. Four of the 2026-09-09 freeze
    (icarm-503, -319, -45, -374) had fewer than 30 good primes below 200 and matched nothing."""
    t = known_curves.load()
    missed = [c["id"] for c in t["curves"] if known_curves.match(c["fingerprint"], t) is None]
    assert missed == [], missed


def _frozen_ainvs(curve_id):
    t = known_curves.load()
    c = next(c for c in t["curves"] if c["id"] == curve_id)
    return "E " + " ".join(c["ainvs"]) + "\n"


@pytest.mark.parametrize("curve_id", ["icarm-503", "icarm-319", "icarm-45", "icarm-374"])
def test_thin_fingerprint_known_curves_are_rejected(curve_id):
    """Known curves with many small bad primes (the self-match failures of the 2026-09 audit)."""
    assert curve_id in expect_known(_frozen_ainvs(curve_id))


def test_post_freeze_icarm_curve_is_rejected():
    """icarm-846 (rank 31, published 2026-09-22) was accepted as new before the 2026-10 refresh."""
    assert "icarm-846" in expect_known(_frozen_ainvs("icarm-846"))


def test_record_curve_is_known():
    msg = expect_known(corpus("rk31.ec"))
    assert "31" in msg


def test_classic_anchor_is_known():
    expect_known(gens_5077a())


def test_baseline_is_novel_and_scores():
    r = verify_text(base5(), FAST_NOVELTY)
    assert r.certified_rank == 5


def test_remodel_of_a_known_curve_is_known():
    # 5077a in the non-reduced model a1 = 2 (change x -> x, y -> y + x): same curve, same a_p.
    out = gpmod.run_gp("E = ellchangecurve(ellinit([0,0,1,-7,6]), [1, 0, 1, 0]);\n"
                       'print("A=", [E.a1, E.a2, E.a3, E.a4, E.a6]);', timeout=30)
    a = [int(x) for x in out["A"].strip("[] ").split(",")]
    expect_known("E " + " ".join(map(str, a)) + "\n")


def test_quadratic_twist_of_a_known_curve_is_known():
    # twist of 5077a by d = -3 via elltwist; |a_p| agree at every good prime.
    out = gpmod.run_gp("E = elltwist(ellinit([0,0,1,-7,6]), -3);\n"
                       'print("A=", [E.a1, E.a2, E.a3, E.a4, E.a6]);', timeout=30)
    a = [int(x) for x in out["A"].strip("[] ").split(",")]
    expect_known("E " + " ".join(map(str, a)) + "\n")


def test_isogenous_curve_of_a_known_curve_is_known():
    # 11a1 has 5-isogenies; ellisomat lists the class as short models [A, B] over Q, which
    # ellminimalmodel turns back into integral minimal models.
    out = gpmod.run_gp("L = ellisomat(ellinit([0,-1,1,-10,-20]), 0, 1)[1];\n"
                       'print("N=", #L); for(i = 1, #L, M = ellminimalmodel(ellinit([0, 0, 0, L[i][1], L[i][2]])); '
                       'print("A", i, "=", [M.a1, M.a2, M.a3, M.a4, M.a6]));', timeout=60)
    n = int(out["N"])
    assert n >= 2
    for i in range(1, n + 1):
        a = [int(x) for x in out[f"A{i}"].strip("[] ").split(",")]
        expect_known("E " + " ".join(map(str, a)) + "\n")


def test_fingerprint_match_requires_enough_good_primes():
    n = len(fp_mod.PRIMES)
    a = [1] * n
    b = [-1] * n
    assert not fp_mod.matches(a, b)
    c = list(a); c[5] = 2
    assert not fp_mod.matches(a, c)
    assert fp_mod.matches(a, a)


def test_candidate_bad_at_most_of_the_known_curves_good_primes_still_matches():
    """A twist by a smooth d is bad at every prime dividing d. Those primes are skipped; the
    match must stand on the primes that remain. The 2026-10-05 patch required compared >= half
    of the known curve's good primes, which is an extra condition and so let such twists
    through (see test_smooth_twist_of_a_known_j1728_curve_is_known)."""
    n = len(fp_mod.PRIMES)
    known = [1] * n
    sparse = [1] * 40 + [-1] * (n - 40)          # 40 compared: over 30, well under n/2
    assert fp_mod.matches(sparse, known)
    too_sparse = [1] * 29 + [-1] * (n - 29)      # 29 compared: under the floor
    assert not fp_mod.matches(too_sparse, known)
    disagree = list(sparse); disagree[7] = 2    # one disagreement at a mutually good prime
    assert not fp_mod.matches(disagree, known)
    assert not hasattr(fp_mod, "MIN_COMPARED_FRACTION")


def _smooth_twist_of(curve_id: str, n_primes: int) -> tuple[str, list[int]]:
    """y^2 = x^3 + a4 d^2 x: the quadratic twist by d = product of the known curve's n_primes
    smallest good primes. val_p(Delta) = 6 at each p | d, so the model is minimal and bad there.
    Returns (text, ainvs)."""
    t = known_curves.load()
    c = next(c for c in t["curves"] if c["id"] == curve_id)
    a = [int(x) for x in c["ainvs"]]
    assert a[0] == a[1] == a[2] == a[4] == 0, "needs a pure y^2 = x^3 + a4 x curve"
    good = [p for p, v in zip(fp_mod.PRIMES, c["fingerprint"]) if v >= 0]
    d = 1
    for p in good[:n_primes]:
        d *= p
    ainvs = [0, 0, 0, a[3] * d * d, 0]
    return (f"E 0 0 0 {ainvs[3]} 0\n", ainvs)


def test_smooth_twist_of_a_known_j1728_curve_is_known():
    """icarm-51 (rank 13, 26-digit a4, 157 good primes below 1000) twisted by the product of its
    83 smallest good primes: a 396-digit a4 (inside the 400-digit cap), minimal at every small
    prime, bad at 83 primes where icarm-51 is good, so 74 primes are compared and all agree.
    The 2026-10-05 rule (compared >= 0.5 * 157) accepted it as new."""
    text, ainvs = _smooth_twist_of("icarm-51", 83)
    assert len(text.split()[4]) <= 400
    msg = expect_known(text)
    assert "rank >= 13" in msg
    fp = fp_mod.fingerprint(ainvs)
    assert fp_mod.good_primes(fp) == 74
    icarm51 = next(c for c in known_curves.load()["curves"] if c["id"] == "icarm-51")
    assert fp_mod.matches(fp, icarm51["fingerprint"])


def test_blind_search_rank6_is_novel_and_certifies():
    """The ladder is real: a curve found by a three-minute blind search passes every gate and
    certifies one rank above the baseline."""
    from util import fixture
    r = verify_text(fixture("rank6_blind.ec"), FAST_NOVELTY)
    assert r.certified_rank == 6


def test_every_frozen_curve_has_enough_good_primes_for_the_candidate_floor():
    """Soundness of the gate against twists by a smooth d (fingerprint.py docstring): a class
    member that passes the candidate floor is good at >= CANDIDATE_MIN_GOOD primes, at most
    (168 - known_good) of which are bad for the known curve, so it is compared at
    >= CANDIDATE_MIN_GOOD - (168 - known_good) primes. That is >= MIN_COMPARED for every frozen
    curve iff known_good >= 168 - CANDIDATE_MIN_GOOD + MIN_COMPARED (= 128). The 2026-10-15
    freeze bottoms out at 139 (icarm-374). A refresh that breaks this must raise CANDIDATE_MIN_GOOD."""
    n = len(fp_mod.PRIMES)
    floor = n - fp_mod.CANDIDATE_MIN_GOOD + fp_mod.MIN_COMPARED
    thin = [(c["id"], fp_mod.good_primes(c["fingerprint"]))
            for c in known_curves.load()["curves"] if fp_mod.good_primes(c["fingerprint"]) < floor]
    assert thin == [], thin
    assert min(fp_mod.good_primes(c["fingerprint"]) for c in known_curves.load()["curves"]) >= 139


def test_candidate_with_too_few_good_primes_is_rejected(monkeypatch):
    """No curve under the 400-digit cap can reach this (a twist's a4 is >= d^2), so the gate
    is exercised with a synthetic fingerprint: the gate must fire before the known-set lookup."""
    n = len(fp_mod.PRIMES)
    synthetic = [1] * (fp_mod.CANDIDATE_MIN_GOOD - 1) + [-1] * (n - fp_mod.CANDIDATE_MIN_GOOD + 1)
    monkeypatch.setattr(fp_mod, "fingerprint", lambda ainvs, gp_bin=None: synthetic)
    with pytest.raises(VerifyError) as ei:
        verify_text(base5(), FAST_NOVELTY)
    assert ei.value.code is ErrorCode.CURVE_TOO_MANY_BAD_PRIMES, ei.value
    assert str(fp_mod.CANDIDATE_MIN_GOOD) in ei.value.message
    ok = [1] * fp_mod.CANDIDATE_MIN_GOOD + [-1] * (n - fp_mod.CANDIDATE_MIN_GOOD)
    monkeypatch.setattr(fp_mod, "fingerprint", lambda ainvs, gp_bin=None: ok)
    assert verify_text(base5(), FAST_NOVELTY).certified_rank == 5     # exactly at the floor: passes


def test_every_frozen_curve_has_a_torsion_key_and_verification():
    for c in known_curves.load()["curves"]:
        assert c["torsion_key"] is not None and c["verification"]["status"] in {"certified", "failed", "unverified"}


def test_bad_prime_floor_holds_with_novelty_off(monkeypatch):
    """--no-novelty skips the known-curve lookup, not the good-prime floor: the floor is what
    makes fingerprint matching (DUPLICATE_CURVE too) sound, and a search loop must not accept a
    curve the official run rejects."""
    from util import FAST
    n = len(fp_mod.PRIMES)
    synthetic = [1] * (fp_mod.CANDIDATE_MIN_GOOD - 1) + [-1] * (n - fp_mod.CANDIDATE_MIN_GOOD + 1)
    monkeypatch.setattr(fp_mod, "fingerprint", lambda ainvs, gp_bin=None: synthetic)
    with pytest.raises(VerifyError) as ei:
        verify_text(base5(), FAST)
    assert ei.value.code is ErrorCode.CURVE_TOO_MANY_BAD_PRIMES, ei.value
