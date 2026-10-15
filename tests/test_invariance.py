"""Facts the audit relied on, frozen as regressions: what cannot move the score."""

import json
import math
import re

import pytest

from conftest import needs_gp
from util import FAST, ROOT, corpus, count_points, gens_5077a

from ec_verify import ErrorCode, VerifyError, verify_text
from ec_verify import gp as gpmod
from ec_verify.parse import parse

pytestmark = needs_gp


def _pt(x, y):
    return gpmod.gp_vec([gpmod.gp_rat(x.numerator, x.denominator), gpmod.gp_rat(y.numerator, y.denominator)])


def _transform_points(ainvs, points, n=None, scale=None):
    """Return P -> n*P (or the image under the u=1/scale change of curve) via gp."""
    vec = gpmod.gp_vec([gpmod.gp_int(a) for a in ainvs])
    body = [f"E = ellinit({vec});"]
    for i, (x, y) in enumerate(points):
        if n is not None:
            body.append(f"Q = ellmul(E, {_pt(x, y)}, {n});")
        else:
            body.append(f"Q = ellchangepoint({_pt(x, y)}, [1/{scale}, 0, 0, 0]);")
        body.append(f'print("X{i}=", numerator(Q[1]), "/", denominator(Q[1])); print("Y{i}=", numerator(Q[2]), "/", denominator(Q[2]));')
    out = gpmod.run_gp("\n".join(body), timeout=120)
    return [(out[f"X{i}"].replace(" ", ""), out[f"Y{i}"].replace(" ", "")) for i in range(len(points))]


def _file(ainvs, pts):
    return "E " + " ".join(str(a) for a in ainvs) + "\n" + "".join(f"P {x} {y}\n" for x, y in pts)


def test_rescaling_above_the_minimality_bound_cannot_move_the_score():
    p = 1000003
    base = verify_text(gens_5077a(), FAST)
    a1, a2, a3, a4, a6 = 0, 0, 1, -7, 6
    scaled = (a1 * p, a2 * p ** 2, a3 * p ** 3, a4 * p ** 4, a6 * p ** 6)
    pts = _transform_points((a1, a2, a3, a4, a6), parse(gens_5077a()).points, scale=p)
    r = verify_text(_file(scaled, pts), FAST)
    assert r.certified_rank == base.certified_rank == 3
    assert r.regulator == pytest.approx(base.regulator, rel=1e-12)
    # H is invariant under u = ±1 only; a u = p re-model multiplies c4^3, c6^2 by p^12, so H grows by log p.
    assert r.naive_height == pytest.approx(base.naive_height + math.log(p), abs=1e-5)
    assert r.naive_height > base.naive_height
    assert r.torsion_order == base.torsion_order


def test_record_witness_survives_multiplying_every_point_by_five():
    text = corpus("rk30.ec")
    s = parse(text)
    pts = _transform_points(s.ainvs, s.points, n=5)
    assert max(len(t) for x, y in pts for t in (x + "/" + y).split("/")) < 2000
    base = verify_text(text, FAST)
    r = verify_text(_file(s.ainvs, pts), FAST)
    assert r.certified_rank == 30
    assert r.regulator == pytest.approx(base.regulator * 5 ** 60, rel=1e-6)
    assert r.height_precision_agreement <= 1e-40


def test_torsion_translate_of_a_generator_is_dependent():
    text = corpus("rkeq10.ec")
    s = parse(text)
    # 2-torsion point of the rkeq10 curve (Dujella's page); frozen here so the test needs no cache.
    T = ("1022753212", "-511376606")
    vec = gpmod.gp_vec([gpmod.gp_int(a) for a in s.ainvs])
    x, y = s.points[0]
    out = gpmod.run_gp(f"E = ellinit({vec});\nT = [{T[0]},{T[1]}];\nprint(\"OT=\", ellorder(E, T));\n"
                       f"Q = elladd(E, {_pt(x, y)}, T);\nprint(\"QX=\", numerator(Q[1]), \"/\", denominator(Q[1]));\n"
                       f"print(\"QY=\", numerator(Q[2]), \"/\", denominator(Q[2]));", timeout=60)
    assert out["OT"] == "2"
    qx, qy = out["QX"].replace(" ", ""), out["QY"].replace(" ", "")
    with pytest.raises(VerifyError) as ei:
        verify_text(text + f"P {qx} {qy}\n", FAST)
    assert ei.value.code is ErrorCode.POINTS_DEPENDENT
    assert "P11" in ei.value.message and "E(Q)_tors" in ei.value.message
    # the torsion point itself is simply dropped
    r = verify_text(text + f"P {T[0]} {T[1]}\n", FAST)
    assert r.certified_rank == 10 and r.points_torsion_dropped == 1 and r.torsion_order == 2


@pytest.mark.parametrize("name", sorted(p.name for p in (ROOT / "tests" / "corpus").glob("*.ec")))
def test_every_corpus_curve_passes_the_height_tripwires(name):
    text = corpus(name)
    r = verify_text(text, FAST)
    assert r.certified_rank == count_points(text)
    assert r.height_precision_agreement <= 1e-40
