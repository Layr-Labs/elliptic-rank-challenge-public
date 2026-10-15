"""Pure formula layer (docs/verifier-architecture.md §4): no gp."""

import pytest

from ec_verify import cells


@pytest.mark.parametrize("structure,key", [([], "1"), ([5], "5"), ([2, 2], "2x2"), ([6, 2], "6x2"), ([8, 2], "8x2")])
def test_torsion_key_follows_pari_order(structure, key):
    assert cells.torsion_key(structure) == key


def test_torsion_key_rejects_non_mazur_shapes():
    with pytest.raises(ValueError):
        cells.torsion_key([2, 2, 2])
    with pytest.raises(ValueError):
        cells.torsion_key([2, 6])          # PARI lists the larger factor first


def test_cell_key():
    assert cells.cell_key("2", 9) == "2/9"
    assert cells.cell_key("1", 0) == "1/0"


def test_contribution_against_a_baseline():
    assert cells.contribution(12.0, 16.0) == pytest.approx(25.0)
    assert cells.contribution(16.0, 16.0) == 0.0            # strict: the record itself scores nothing
    assert cells.contribution(20.0, 16.0) == 0.0


def test_breakthrough_is_above_100_and_strictly_decreasing_in_H():
    assert cells.contribution(0.0, None) == 200.0
    assert cells.contribution(cells.H_REF, None) == pytest.approx(150.0)
    assert cells.contribution(100.0, None) == pytest.approx(100 + 100 / 3)
    a, b = cells.contribution(40.0, None), cells.contribution(40.000001, None)
    assert 100.0 < b < a


def test_breakthrough_keeps_its_slope_for_huge_curves():
    # The 400-digit cap is on the raw a-invariants: a 400-digit a1 gives H above 900. There is no
    # clamp, so a smaller curve always scores more, and nothing ever falls to 100.
    assert 100.0 < cells.contribution(1e6, None) < cells.contribution(1e3, None) < cells.contribution(460.0, None)


def test_contribution_is_strictly_decreasing_everywhere_it_is_positive():
    for base in (None, 30.0):
        hs = [1.0, 5.0, 10.0, 20.0, 29.999999]
        vals = [cells.contribution(h, base) for h in hs]
        assert vals == sorted(vals, reverse=True) and len(set(vals)) == len(vals)


def test_board_total_sums_in_sorted_key_order():
    c = {"2/9": 0.1, "1/11": 0.2, "6x2/4": 0.3}
    assert cells.board_total(c) == (0.2 + 0.1) + 0.3          # "1/11", "2/9", "6x2/4"
    assert cells.board_total({}) == 0.0


def test_cell_rank_and_the_scored_rank_floor():
    assert cells.cell_rank("6x2/4") == 4
    assert cells.cell_rank("1/0") == 0
    assert cells.MIN_SCORED_RANK == 1


def test_curve_contribution_is_zero_at_rank_zero_and_the_cell_formula_otherwise():
    assert cells.curve_contribution(0, 1.0, 20.0) == 0.0 and cells.curve_contribution(0, 1.0, None) == 0.0
    assert cells.curve_contribution(5, 12.0, 16.0) == cells.contribution(12.0, 16.0)
    assert cells.curve_contribution(1, 40.0, None) == cells.contribution(40.0, None)
