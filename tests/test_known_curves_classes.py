"""known_curves.classes / representative: the frozen set grouped into isogeny classes per rank.

Isogenous curves share conductor and L-function, so their |a_p| fingerprints coincide and the
novelty gate treats them as one known curve; anything that counts "known curves of rank r"
(the per-torsion record table) must do the same, or a class
uploaded whole counts once per member.
"""

from ec_verify import known_curves

FP_A = [1, 2, -1, 3] + [1] * 30  # pad to reach MIN_COMPARED threshold
FP_B = [1, 2, 0, 3] + [2] * 30


def _known(*curves):
    return {"version": 2, "frozen_at": "2026-10-05T05:49:16+00:00", "curves": list(curves)}


def _c(cid, rank, ainvs, fp):
    return {"id": cid, "source": cid.split("-")[0], "source_id": cid.split("-", 1)[1], "rank_claimed": rank,
            "ainvs": [str(a) for a in ainvs], "fingerprint": fp, "duplicate_of_earlier": False,
            "torsion": None, "torsion_key": None, "naive_height": 0.0,
            "verification": {"status": "unverified", "certified_rank": None, "code": None,
                             "verified_at": None, "pari_version": None, "verifier_version": None}}


def test_isogenous_curves_form_one_class_with_the_smallest_member_as_representative():
    big = _c("icarm-2", 3, [0, 0, 0, -100000, 1000000], FP_A)
    small = _c("icarm-1", 3, [0, 0, 1, -7, 6], FP_A)
    classes = known_curves.classes(_known(big, small))
    assert list(classes) == [(tuple(FP_A), 3)]
    members = classes[(tuple(FP_A), 3)]
    assert [m["id"] for m in members] == ["icarm-2", "icarm-1"]          # input order kept
    assert known_curves.representative(members)["id"] == "icarm-1"


def test_same_fingerprint_different_rank_stays_separate():
    a = _c("icarm-1", 3, [0, 0, 1, -7, 6], FP_A)
    tw = _c("icarm-9", 12, [0, 0, 0, -7 * 9, 6 * 27], FP_A)          # twist at another rank
    classes = known_curves.classes(_known(a, tw))
    assert sorted(classes) == [(tuple(FP_A), 3), (tuple(FP_A), 12)]


def test_different_fingerprint_same_rank_stays_separate():
    a = _c("icarm-5", 3, [0, 0, 1, -7, 6], FP_A)
    b = _c("icarm-4", 3, [0, 0, 0, -1, 1], FP_B)
    classes = known_curves.classes(_known(a, b))
    assert len(classes) == 2
    assert all(len(m) == 1 for m in classes.values())


def test_tie_on_height_picks_the_smaller_id_numerically():
    a = _c("icarm-10", 3, [0, 0, 1, -7, 6], FP_A)
    b = _c("icarm-9", 3, [0, 0, 1, -7, 6], FP_A)
    assert known_curves.representative([a, b])["id"] == "icarm-9"       # 9 < 10, not "icarm-10" < "icarm-9"


def test_the_frozen_set_has_5586_curves_in_4239_classes_covering_every_record_once():
    t = known_curves.load()
    classes = known_curves.classes(t)
    assert len(t["curves"]) == 5586
    assert len(classes) == 4239
    members = [m["id"] for ms in classes.values() for m in ms]
    assert len(members) == len(set(members)) == 5586
    # the gate agrees with the grouping: members of one class match each other's fingerprint
    for ms in list(classes.values())[:50]:
        for m in ms[1:]:
            assert known_curves.match(m["fingerprint"], {"curves": [ms[0]]}) is ms[0]


def test_validate_ledger_requires_v2_fields():
    bad = _c("icarm-1", 3, [0, 0, 1, -7, 6], FP_A)
    bad.pop("verification")
    import pytest
    with pytest.raises(ValueError, match="verification"):
        known_curves.validate_ledger(_known(bad))


def test_validate_ledger_reports_a_malformed_verification_block_as_a_value_error():
    # A missing certified_rank or status must name the entry in a ValueError, not surface as a bare
    # KeyError at the end of a full rebuild.
    import pytest
    bad = _c("icarm-1", 3, [0, 0, 1, -7, 6], FP_A)
    bad["verification"] = {"status": "certified"}
    with pytest.raises(ValueError, match="icarm-1"):
        known_curves.validate_ledger(_known(bad))
    bad["verification"] = {"certified_rank": 3}
    with pytest.raises(ValueError, match="icarm-1"):
        known_curves.validate_ledger(_known(bad))
