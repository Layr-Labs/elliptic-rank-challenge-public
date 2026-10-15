# tests/test_records.py
"""ec_verify/records.py: the frozen cell table (docs/records.md §1.3) and its invariants."""

import json
import pathlib

import pytest

from ec_verify import records

FIX = pathlib.Path(__file__).resolve().parent / "fixtures" / "records_fixture.json"


def test_fixture_loads_and_validates():
    t = records.load(FIX)
    records.validate(t)
    assert t["version"] == 1 and t["scoring_version"] == 6


def test_lookups():
    t = records.load(FIX)
    assert records.baseline(t, "1/5") == 4.0
    assert records.record(t, "1/5")["id"] == "fixture-r5"
    assert records.baseline(t, "2/10") == 12.0            # the baseline is the record height, frozen
    assert records.record(t, "2/10")["naive_height"] == 12.0
    assert records.cell(t, "3/7") is None
    assert records.baseline(t, "3/7") is None and records.record(t, "3/7") is None   # a breakthrough cell


@pytest.mark.parametrize("mutate,msg", [
    (lambda t: t["cells"]["1/5"].pop("record"), "record"),
    (lambda t: t["cells"].__setitem__("bogus", t["cells"]["1/5"]), "key"),
    (lambda t: t["cells"]["1/5"]["record"].__setitem__("naive_height", -1.0), "positive"),
    (lambda t: t["cells"]["1/5"]["record"].__setitem__("naive_height", 0.0), "positive"),    # contribution divides by it
    (lambda t: t.__setitem__("scoring_version", 5), "scoring_version"),       # a v5 table (baselines) is refused
    (lambda t: t["cells"]["1/5"]["record"].__setitem__("certified_rank", 4), "certified_rank"),   # cannot reach rank 5
])
def test_validate_rejects_malformed_tables(mutate, msg):
    t = json.loads(FIX.read_text())
    mutate(t)
    with pytest.raises(ValueError, match=msg):
        records.validate(t)


def test_shipped_table_validates_and_matches_the_scoring_version():
    t = records.load()
    records.validate(t)
    from ec_verify.score import SCORING_VERSION
    assert t["scoring_version"] == SCORING_VERSION


def test_frozen_table_covers_every_torsion_group_and_matches_the_ledger():
    from ec_verify import cells, known_curves
    t = records.load()
    led = known_curves.load()
    assert t["version"] >= 1 and t["ledger_frozen_at"] == led["frozen_at"]
    certified = [c for c in led["curves"] if c["verification"]["status"] == "certified"]
    assert len(certified) >= 3000
    # A curve certified at rank R reaches every rank 1..R (list fewer of its points), so the record
    # of cell (T, r) is the smallest certified curve of torsion T with certified rank >= r.
    expected = {}
    for c in certified:
        R = int(c["verification"]["certified_rank"])
        for r in range(cells.MIN_SCORED_RANK, R + 1):      # rank 0 proves nothing: never a cell
            k = cells.cell_key(c["torsion_key"], r)
            expected[k] = min(expected.get(k, 1e9), c["naive_height"])
    assert set(t["cells"]) == set(expected)
    for k, c in t["cells"].items():
        assert c["record"]["naive_height"] == expected[k]
        assert c["record"]["certified_rank"] >= c["rank"] and c["exact_count"] <= c["known_count"]
        assert "filled_from_rank" not in c
    assert {c["torsion_key"] for c in t["cells"].values()} == {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "12", "2x2", "4x2", "6x2", "8x2"}
    assert "h_cap" not in t and not any("baseline" in c for c in t["cells"].values())


def test_frozen_table_has_no_rank_gap_below_a_torsion_groups_top_rank():
    # The verifier awards a breakthrough to any key absent from the table; the documented rule is
    # "no known curve of that torsion reaches the rank", so every rank from 1 up to the group's
    # maximum must be present (build_records fills gaps from the next certified rank up).
    from ec_verify import cells
    t = records.load()
    by_group = {}
    for c in t["cells"].values():
        by_group.setdefault(c["torsion_key"], set()).add(int(c["rank"]))
    for tkey, ranks in by_group.items():
        assert ranks == set(range(cells.MIN_SCORED_RANK, max(ranks) + 1)), tkey


def test_icarm_height_is_twelve_times_ours():
    import json, pathlib
    raw = json.loads((pathlib.Path(__file__).resolve().parent.parent / "tools" / "known_raw.json").read_text())
    from ec_verify.curve import naive_height
    worst = max(abs(c["naive_height_source"] / 12 - naive_height(c["ainvs"]))
                for c in raw["curves"] if c.get("naive_height_source") is not None)
    assert worst < 1e-5
