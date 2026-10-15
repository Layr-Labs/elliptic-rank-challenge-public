"""gp-backed gates: minimality at small primes, torsion stripping, claims."""

import pytest

from conftest import needs_gp
from util import FAST, gens_5077a

from ec_verify import ErrorCode, VerifyError, verify_text

pytestmark = needs_gp


def expect(text, code):
    with pytest.raises(VerifyError) as ei:
        verify_text(text, FAST)
    assert ei.value.code is code, ei.value


def test_non_minimal_model_rejected():
    # 5077a scaled by u = 1/5: [0,0,125,-4375,93750] with (-2,3) -> (-50,375); Δ picks up 5^12.
    expect("E 0 0 125 -4375 93750\nP -50 375\n", ErrorCode.CURVE_NON_MINIMAL)


def test_torsion_only_certifies_zero():
    # 11a1: (5,5) has order 5. It is on the curve but never counts.
    r = verify_text("E 0 -1 1 -10 -20\nP 5 5\n", FAST)
    assert r.certified_rank == 0
    assert r.torsion_order == 5
    assert r.points_torsion_dropped == 1
    assert r.torsion == [5] and r.torsion_key == "5" and r.cell == "5/0"


def test_torsion_mixed_with_free_points():
    # 37a1 (rank 1, trivial torsion) generator (0,0); add nothing torsion -> 1. Then 11a1 torsion + free? 11a1 has rank 0.
    r = verify_text("E 0 0 1 -1 0\nP 0 0\n", FAST)
    assert r.certified_rank == 1 and r.torsion_order == 1
    assert r.cell == "1/1"


def test_claim_is_recorded_not_trusted():
    r = verify_text(gens_5077a().replace("#CLAIM rank=3", "#CLAIM rank=7"), FAST)
    assert r.certified_rank == 3
    assert r.claimed_rank == 7 and r.claim_discrepancy
    assert "not established" in r.verified_claim


def test_no_points_scores_zero_rank():
    r = verify_text("E 0 0 1 -7 6\n", FAST)
    assert r.certified_rank == 0 and r.regulator == 1.0


def test_gp_version_mismatch_code(monkeypatch):
    from ec_verify import gp as gpmod
    monkeypatch.setattr(gpmod, "PINNED_VERSION", "0.0.0")
    expect(gens_5077a(), ErrorCode.GP_VERSION_MISMATCH)


def test_strip_torsion_reports_the_structure():
    from ec_verify import torsion
    order, structure, keep, dropped = torsion.strip_torsion((0, -1, 1, -10, -20), [], gp_bin=None)
    assert (order, structure) == (5, [5])
    order, structure, keep, dropped = torsion.strip_torsion((0, 0, 0, -1, 0), [], gp_bin=None)
    assert (order, structure) == (4, [2, 2])
    order, structure, keep, dropped = torsion.strip_torsion((1, 1, 1, -5, 2), [], gp_bin=None)
    assert (order, structure) == (8, [4, 2])
    order, structure, keep, dropped = torsion.strip_torsion((0, 0, 1, -7, 6), [], gp_bin=None)
    assert (order, structure) == (1, [])


def test_parse_gp_ints():
    import pytest
    from ec_verify.gp import GPError
    from ec_verify.torsion import parse_gp_ints
    assert parse_gp_ints("[6, 2]") == [6, 2]
    assert parse_gp_ints("[]") == []
    assert parse_gp_ints("[5]") == [5]
    with pytest.raises(GPError):
        parse_gp_ints("[6, x]")
    for malformed in ("[6, ]", "[6,,2]", "[,]"):           # an empty entry is not an integer: fail closed
        with pytest.raises(GPError):
            parse_gp_ints(malformed)


def test_result_carries_cell_and_sizes():
    from util import base5
    r = verify_text(base5(), FAST)
    assert r.cell == "1/5" and r.torsion == [] and r.torsion_key == "1"
    assert r.naive_height == pytest.approx(3.422446) and r.naive_height_source == pytest.approx(12 * 3.422446, abs=1e-5)
    assert len(r.fingerprint) == 168


def test_unrecognised_torsion_structure_from_gp_is_a_gp_error_not_a_traceback(monkeypatch):
    # cells.torsion_key raises ValueError for a shape outside Mazur's list; from gp that is a
    # malformed STRUCT line and must surface as the documented GP_ERROR code.
    import ec_verify.verify as vmod
    real = vmod.torsion_mod.strip_torsion

    def odd(ainvs, pts, gp_bin=None):
        tors, structure, free, dropped = real(ainvs, pts, gp_bin=gp_bin)
        return tors, [2, 6], free, dropped

    monkeypatch.setattr(vmod.torsion_mod, "strip_torsion", odd)
    with pytest.raises(VerifyError) as ei:
        verify_text(gens_5077a(), FAST)
    assert ei.value.code is ErrorCode.GP_ERROR and "[2, 6]" in ei.value.message
