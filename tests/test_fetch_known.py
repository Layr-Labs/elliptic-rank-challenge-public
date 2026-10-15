"""tools/fetch_known.py: raw mirror v2, strict parsing, idempotent merge."""

import json
import pathlib
import sys

import pytest

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "tools"))
import fetch_known  # noqa: E402
import fetch_lmfdb  # noqa: E402
import fetch_records  # noqa: E402

GOOD = {"id": 7, "ainvs": ["0", "0", "1", "-7", "6"], "rank_lower_bound": 3, "torsion": [],
        "points": [["1", "2"], ["3/4", "5/8"], ["5", "6"]], "naive_height": 17.45,
        "created_at": "2026-01-01 00:00:00", "updated_at": "2026-02-01 00:00:00"}


def test_good_record_parses_to_v2_row():
    r = fetch_known.icarm_row(GOOD)
    assert r["id"] == "icarm-7" and r["source_id"] == "7" and r["ainvs"] == [0, 0, 1, -7, 6]
    assert r["points"] == [["1", "2"], ["3/4", "5/8"], ["5", "6"]] and r["rank_claimed"] == 3
    assert r["torsion_source"] == [] and r["naive_height_source"] == 17.45
    assert r["source_updated_at"] == "2026-02-01 00:00:00"


def test_bracketed_string_ainvs_parse_too():
    r = fetch_known.icarm_row({**GOOD, "ainvs": "[0, 0, 1, -7, 6]"})
    assert r["ainvs"] == [0, 0, 1, -7, 6]


@pytest.mark.parametrize("bad", [
    {**GOOD, "ainvs": "[0, 0, 1, -7]"}, {**GOOD, "ainvs": "[0, 0, 1, -7/2, 6]"}, {k: v for k, v in GOOD.items() if k != "id"},
    {**GOOD, "points": [["1"]]}, {**GOOD, "torsion": [2, 2, 2]},
])
def test_malformed_record_raises(bad):
    with pytest.raises(fetch_known.BadRecord):
        fetch_known.icarm_row(bad)


@pytest.mark.parametrize("ainvs", [[0, 0, 0, -2.5, 1], [0, 0, 0, True, 1], ["0", "0", "0", "-2.5", "1"], [0, 0, 0, 2.0, 1]])
def test_non_integral_ainvs_are_rejected_not_truncated(ainvs):
    # int(-2.5) == -2 would mirror a different curve than the published one, and the published
    # one would then pass the KNOWN_CURVE gate as "new".
    with pytest.raises(fetch_known.BadRecord, match="integers"):
        fetch_known.icarm_row({**GOOD, "ainvs": ainvs})


def test_icarm_rows_fails_loudly_instead_of_dropping(capsys):
    with pytest.raises(SystemExit):
        fetch_known.icarm_rows([GOOD, {**GOOD, "id": 8, "ainvs": "[1]"}])
    assert "id 8" in capsys.readouterr().err


def test_merge_is_idempotent_and_keeps_verification_for_unchanged_rows():
    row = fetch_known.icarm_row(GOOD)
    raw = fetch_known.merge({"version": 2, "curves": []}, [row], now="2026-10-08T00:00:00+00:00")
    assert raw["curves"][0]["first_seen"] == "2026-10-08T00:00:00+00:00"
    raw["curves"][0]["verification"] = {"status": "certified", "certified_rank": 3}
    again = fetch_known.merge(raw, [row], now="2026-10-09T00:00:00+00:00")
    assert again["curves"][0]["verification"]["status"] == "certified"
    assert again["curves"][0]["first_seen"] == "2026-10-08T00:00:00+00:00"
    updated = fetch_known.merge(again, [fetch_known.icarm_row({**GOOD, "updated_at": "2026-03-01 00:00:00"})], now="2026-10-10T00:00:00+00:00")
    assert "verification" not in updated["curves"][0] and updated["curves"][0]["first_seen"] == "2026-10-08T00:00:00+00:00"


def test_merge_replaces_an_entry_whose_content_changed_without_a_timestamp():
    # Dujella/LMFDB rows have no updated_at: a refreshed page with an extra point (or a corrected
    # rank/torsion) must replace the entry and drop its stale verification, keeping first_seen.
    row = {"id": "dujella-rk20", "source": "dujella", "source_id": "rk20", "ainvs": [1, 0, 0, -1, 1], "points": [["1", "1"]],
           "rank_claimed": 20, "torsion_source": None, "naive_height_source": None, "source_created_at": None, "source_updated_at": None}
    raw = fetch_known.merge({"version": 2, "curves": []}, [row], now="t1")
    raw["curves"][0]["verification"] = {"status": "certified", "certified_rank": 1}
    same = fetch_known.merge(raw, [dict(row)], now="t2")
    assert same["curves"][0]["verification"]["status"] == "certified"
    more = fetch_known.merge(raw, [{**row, "points": [["1", "1"], ["2", "3"]]}], now="t3")
    assert "verification" not in more["curves"][0] and more["curves"][0]["points"] == [["1", "1"], ["2", "3"]]
    assert more["curves"][0]["first_seen"] == "t1"
    retors = fetch_known.merge(raw, [{**row, "torsion_source": [2]}], now="t4")
    assert "verification" not in retors["curves"][0]


def test_merge_sorts_by_source_then_numeric_id():
    rows = [fetch_known.icarm_row({**GOOD, "id": 10}), fetch_known.icarm_row({**GOOD, "id": 9}),
            {"id": "dujella-rk20", "source": "dujella", "source_id": "rk20", "ainvs": [1, 0, 0, -1, 1], "points": [],
             "rank_claimed": 20, "torsion_source": None, "naive_height_source": None, "source_created_at": None, "source_updated_at": None}]
    raw = fetch_known.merge({"version": 2, "curves": []}, rows, now="t")
    assert [c["id"] for c in raw["curves"]] == ["icarm-9", "icarm-10", "dujella-rk20"]


def test_merge_keeps_entries_that_vanished_upstream():
    row = fetch_known.icarm_row(GOOD)
    raw = fetch_known.merge({"version": 2, "curves": []}, [row], now="t1")
    other = fetch_known.icarm_row({**GOOD, "id": 8})
    again = fetch_known.merge(raw, [other], now="t2")        # id 7 no longer in the fetch
    assert [c["id"] for c in again["curves"]] == ["icarm-7", "icarm-8"]
    assert again["curves"][0]["first_seen"] == "t1"


def test_merge_fills_first_seen_for_a_changed_entry_that_lacked_it():
    # A hand-edited or v1-era mirror entry without first_seen must not abort the whole fetch with a
    # KeyError when its upstream content changes; it gets `now`, like a new entry.
    row = fetch_known.icarm_row(GOOD)
    raw = fetch_known.merge({"version": 2, "curves": []}, [row], now="t1")
    raw["curves"][0].pop("first_seen")
    changed = fetch_known.merge(raw, [fetch_known.icarm_row({**GOOD, "updated_at": "2026-03-01 00:00:00"})], now="t2")
    assert changed["curves"][0]["first_seen"] == "t2"


# ---- Dujella per-torsion pages: several curves per page, torsion points listed, sign typos ----

TORS_PAGE = """<html><body><h1>Torsion group Z/2Z &times; Z/2Z, rank = 2</h1>
<p>Someone (2001)</p>
<p>y<sup>2</sup> + xy = x<sup>3</sup> - x<sup>2</sup> - 4x + 4</p>
<p>Torsion points: O, [1, -1], [2, -2], [-2, 2]</p>
<p>Independent points of infinite order:<br>
P<sub>1</sub> = [3, 2]<br>P<sub>2</sub> = [5,
7]</p>
<p>Someone else (1999)</p>
<p>y<sup>2</sup> + xy + y = x<sup>3</sup> + x<sup>2</sup> - -5x + 3</p>
<p>Torsion points: O, [-3, 1]</p>
<p>Independent points of infinite order: P<sub>1</sub> = [-3, 0]&nbsp; P<sub>2</sub> = [3/4, 1/8]</p>
</body></html>"""


def test_on_curve_uses_exact_rationals():
    assert fetch_records.on_curve([0, 0, 1, -7, 6], "0", "2")
    assert fetch_records.on_curve([0, 0, 0, 18372, -33031], "764/25", "93287/125")
    assert not fetch_records.on_curve([0, 0, 1, -7, 6], "0", "3")


def test_parse_curves_reads_every_curve_on_a_torsion_page_and_skips_torsion_points():
    curves = fetch_records.parse_curves(TORS_PAGE)
    assert [c["ainvs"] for c in curves] == [[1, -1, 0, -4, 4], [1, 1, 1, -5, 3]]
    assert curves[0]["points"] == [("3", "2"), ("5", "7")]
    assert curves[1]["points"] == [("-3", "0"), ("3/4", "1/8")]
    for c in curves:
        assert all(fetch_records.on_curve(c["ainvs"], x, y) for x, y in c["points"])


def test_parse_curves_resolves_a_doubled_sign_by_the_listed_points():
    # "- -5x" on the second curve: only a4 = -5 puts [-3, 0] on the curve.
    assert fetch_records.parse_curves(TORS_PAGE)[1]["ainvs"][3] == -5
    bad = TORS_PAGE.replace("P<sub>1</sub> = [-3, 0]", "P<sub>1</sub> = [-3, 7]")
    with pytest.raises(ValueError, match="not on the curve"):
        fetch_records.parse_curves(bad)


def test_parse_curves_fails_loudly_when_an_equation_has_no_point_list():
    with pytest.raises(ValueError, match="equations on the page"):
        fetch_records.parse_curves(TORS_PAGE + "<p>y<sup>2</sup> = x<sup>3</sup> - 2</p>")


def test_parse_curves_fails_loudly_when_a_point_does_not_parse():
    # a point that is not a pair would silently drop a point and understate the rank
    with pytest.raises(ValueError, match="2 points announced, 1 parsed"):
        fetch_records.parse_curves(TORS_PAGE.replace("[3, 2]", "[3; 2]"))
    # a number broken by markup is re-joined and then fails the on-curve check instead
    with pytest.raises(ValueError, match="not on the curve"):
        fetch_records.parse_curves(TORS_PAGE.replace("[3, 2]", "[3, 2<br>3]"))


def test_parse_curves_reports_problems_as_dicts_when_asked():
    problems = []
    curves = fetch_records.parse_curves(TORS_PAGE.replace("P<sub>1</sub> = [-3, 0]", "P<sub>1</sub> = [-3, 7]"), problems)
    assert len(curves) == 2 and curves[1]["points"] == [("3/4", "1/8")] and curves[1]["listed"] == 2
    assert [p["kind"] for p in problems] == ["dropped"] and problems[0]["index"] == 2
    problems = []
    curves = fetch_records.parse_curves(TORS_PAGE.replace("- -5x", "+ 5x"), problems)   # neither point fits a4 = +5
    assert len(curves) == 1 and problems[0]["kind"] == "omitted" and problems[0]["ainvs"] == [1, 1, 1, 5, 3]


def test_dujella_rows_cover_the_torsion_pages_with_stable_hash_ids(tmp_path, monkeypatch):
    (tmp_path / "z2z2.html").write_text(TORS_PAGE)
    (tmp_path / "rkeq10.html").write_text(TORS_PAGE.replace("rank = 2", "rank = 10"))
    monkeypatch.setattr(fetch_known, "CACHE", tmp_path)
    rows = fetch_known.dujella_rows()
    ids = [r["id"] for r in rows]
    assert ids[0] == "dujella-rkeq10"                      # legacy row: first curve, rank from the file name
    assert rows[0]["rank_claimed"] == 10 and len(rows[0]["points"]) == 2
    assert len(rows) == 2 and ids[1].startswith("dujella-rkeq10-")   # second curve, once, by a-invariant hash
    assert rows[1]["ainvs"] == [1, 1, 1, -5, 3] and rows[1]["rank_claimed"] == 2
    assert rows[1]["points"] == [["-3", "0"], ["3/4", "1/8"]] and rows[1]["torsion_source"] is None
    assert fetch_known.dujella_hash_id("z2z2", [1, 1, 1, -5, 3]) == ids[1].replace("rkeq10", "z2z2")


def test_dujella_rows_prefer_the_historical_id_over_a_hash_id_from_an_earlier_page(tmp_path, monkeypatch):
    # rk20.html sorts before rk24.html; if the rank-24 record is also printed on rk20, it must be
    # mirrored once, as dujella-rk24, not as dujella-rk20-<hash> plus dujella-rk24.
    (tmp_path / "rk20.html").write_text(TORS_PAGE)
    second = TORS_PAGE[TORS_PAGE.index("<p>Someone else"):]
    (tmp_path / "rk24.html").write_text("<html><body>" + second)
    monkeypatch.setattr(fetch_known, "CACHE", tmp_path)
    ids = [r["id"] for r in fetch_known.dujella_rows()]
    assert ids == ["dujella-rk20", "dujella-rk24"]


def test_dujella_rows_abort_on_an_unlisted_misprinted_equation(tmp_path, monkeypatch, capsys):
    (tmp_path / "z5.html").write_text(TORS_PAGE.replace("- -5x", "+ 5x"))
    monkeypatch.setattr(fetch_known, "CACHE", tmp_path)
    with pytest.raises(SystemExit):
        fetch_known.dujella_rows()
    assert "not in DUJELLA_MISPRINTS" in capsys.readouterr().err
    monkeypatch.setattr(fetch_known, "DUJELLA_MISPRINTS", {(1, 1, 1, 5, 3)})
    assert [r["ainvs"] for r in fetch_known.dujella_rows()] == [[1, -1, 0, -4, 4]]


# ---- Cremona ecdata allgens: the LMFDB's curves with their generators ----

ALLGENS = """37 a 1 [0,0,1,-1,0] 1 [] [0:0:1]
43 a 1 [0,1,1,0,0] 1 [] [0:0:1]
11 a 1 [0,-1,1,-10,-20] 0 [5] [5:5:1]
117 a 2 [1,-1,1,-41,96] 1 [2,2] [6:3:1] [-7:3:1] [5:-3:1]
1088 j 1 [0,1,0,-25,39] 2 [2] [1:4:1] [5:8:1] [3:0:1]
3990 r 1 [1,1,1,-8,8] 1 [2] [2:-7:2] [1:-1:1]
"""


def test_allgens_line_parses_label_rank_torsion_and_infinite_order_generators_only():
    c = fetch_lmfdb.parse_line("117 a 2 [1,-1,1,-41,96] 1 [2,2] [6:3:1] [-7:3:1] [5:-3:1]")
    assert c["label"] == "117a2" and c["ainvs"] == [1, -1, 1, -41, 96] and c["rank"] == 1
    assert c["torsion"] == [2, 2] and c["gens"] == [["6", "3"]]
    frac = fetch_lmfdb.parse_line("3990 r 1 [1,1,1,-8,8] 1 [2] [2:-7:2] [1:-1:1]")
    assert frac["gens"] == [["1", "-7/2"]]


def test_select_keeps_the_smallest_curves_per_torsion_and_rank_and_no_rank_zero():
    sel = fetch_lmfdb.select(ALLGENS.splitlines(), per_cell=1)
    assert [c["label"] for c in sel] == ["43a1", "117a2", "1088j1", "3990r1"]      # 43a1 is smaller than 37a1; conductor order
    assert all(c["rank"] >= 1 for c in sel)
    two = fetch_lmfdb.select(ALLGENS.splitlines(), per_cell=2)
    assert {c["label"] for c in two} == {"1088j1", "117a2", "3990r1", "43a1", "37a1"}


def test_lmfdb_rows_carry_points_and_torsion_and_override_the_pointless_anchor(tmp_path, monkeypatch):
    (tmp_path / "anchor_curves.json").write_text(json.dumps([
        {"id": "37a1", "rank": 1, "ainvs": [0, 0, 1, -1, 0], "source": "lmfdb"},
        {"id": "11a1", "rank": 0, "ainvs": [0, -1, 1, -10, -20], "source": "lmfdb"}]))
    (tmp_path / "lmfdb_curves.json").write_text(json.dumps({"curves": [
        {"label": "37a1", "ainvs": [0, 0, 1, -1, 0], "rank": 1, "torsion": [], "gens": [["0", "0"]], "conductor": 37}]}))
    monkeypatch.setattr(fetch_known, "TOOLS", tmp_path)
    rows = {r["id"]: r for r in fetch_known.anchor_rows()}
    assert set(rows) == {"lmfdb-37a1", "lmfdb-11a1"}
    assert rows["lmfdb-37a1"]["points"] == [["0", "0"]] and rows["lmfdb-37a1"]["torsion_source"] == []
    assert rows["lmfdb-11a1"]["points"] == [] and rows["lmfdb-11a1"]["torsion_source"] is None
