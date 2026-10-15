"""Board layer (docs/verifier-architecture.md §4, §5) against the fixture cell table; novelty off (FAST) so the
corpus curves can be used as board members."""

import json
import pathlib

import pytest

from conftest import needs_gp
from util import FAST, FAST_NOVELTY, base5, corpus, fixture, gens_5077a

from ec_verify import ErrorCode, VerifyError, cells
from ec_verify.board import verify_board
from ec_verify.verify import VerifyConfig

pytestmark = needs_gp
TABLE = json.loads((pathlib.Path(__file__).resolve().parent / "fixtures" / "records_fixture.json").read_text())
R6 = fixture("rank6_blind.ec")       # trivial torsion, rank 6, H 4.288426: no fixture cell -> breakthrough
B5 = base5()                         # trivial, rank 5, H 3.422446: fixture baseline 4.0 -> 14.44%
G3 = gens_5077a()                    # trivial, rank 3, H 1.454278: it is the fixture record itself -> 0
Z2R10 = corpus("rkeq10.ec")          # Z/2, rank 10, H 11.510206: fixture record 12.0 -> 4.08%


def test_empty_board_scores_zero_and_lists_known_cells():
    b = verify_board({}, None, FAST, TABLE)
    assert b.board_total == 0.0 and b.curves == [] and b.parent_gate == "skipped"
    assert set(b.cells) == {"1/3", "1/5", "2/10"} and all(c.contribution == 0.0 for c in b.cells.values())


def test_contributions_and_total():
    b = verify_board({"b5.ec": B5, "r6.ec": R6, "z2.ec": Z2R10}, None, FAST, TABLE)
    c15, c16, c210 = b.cells["1/5"], b.cells["1/6"], b.cells["2/10"]
    assert c15.contribution == pytest.approx(100 * (4.0 - 3.422446) / 4.0)
    assert c16.contribution == pytest.approx(cells.contribution(4.288426, None)) and c16.known_record is None
    assert c210.contribution == pytest.approx(100 * (12.0 - 11.510206) / 12.0) and c210.known_record["id"] == "fixture-z2-r10"
    assert b.board_total == pytest.approx(c15.contribution + c16.contribution + c210.contribution)
    js = b.to_json()
    by_file = {c["file"]: c for c in js["curves"]}
    assert by_file["b5.ec"]["holds_cell"] and by_file["b5.ec"]["contribution"] == pytest.approx(c15.contribution)
    assert by_file["z2.ec"]["contribution"] == pytest.approx(100 * (12.0 - 11.510206) / 12.0)
    assert by_file["r6.ec"]["cell_record"] is None and by_file["r6.ec"]["contribution"] == pytest.approx(cells.contribution(4.288426, None))
    assert js["parent_total"] is None and js["improvement"] is None       # no parent: nothing to compare with
    assert {i["cell"] for i in js["improvements"]} == {"1/5", "1/6", "2/10"}


def test_two_curves_in_one_cell_both_listed_smaller_counts(tmp_path):
    # Same cell 1/5: base5 (H 3.42, holds) and rank6_blind with one point removed (rank 5,
    # H 4.288426 > baseline 4.0: holds nothing). Both verify; both are listed; the smaller counts.
    r6_minus = "\n".join(l for l in R6.splitlines() if not l.startswith("P -633 10822")) + "\n"
    b = verify_board({"a.ec": B5, "b.ec": r6_minus}, None, FAST, TABLE)
    assert [f for f, _ in b.curves] == ["a.ec", "b.ec"]
    assert b.cells["1/5"].board_best_file == "a.ec"
    js = {c["file"]: c for c in b.to_json()["curves"]}
    assert js["a.ec"]["holds_cell"] and not js["b.ec"]["holds_cell"]


def test_duplicate_fingerprint_rejected():
    with pytest.raises(VerifyError) as ei:
        verify_board({"a.ec": B5, "b.ec": B5}, None, FAST, TABLE)
    assert ei.value.code is ErrorCode.DUPLICATE_CURVE and "a.ec" in ei.value.message and "b.ec" in ei.value.message


def test_per_file_failure_names_the_file():
    with pytest.raises(VerifyError) as ei:
        verify_board({"ok.ec": B5, "bad.ec": "E 0 0 0 0 0\n"}, None, FAST, TABLE)
    assert ei.value.code is ErrorCode.CURVE_SINGULAR and ei.value.message.startswith("board/bad.ec: ")


def test_too_many_curves():
    files = {f"c{i:03d}.ec": B5 for i in range(cells.MAX_BOARD_FILES + 1)}
    with pytest.raises(VerifyError) as ei:
        verify_board(files, None, FAST, TABLE)
    assert ei.value.code is ErrorCode.TOO_MANY_CURVES


def test_regression_gate_cases():
    parent = {"b5.ec": B5, "r6.ec": R6}
    # unchanged: ok, checked
    b = verify_board(dict(parent), parent, FAST, TABLE)
    assert b.parent_gate == "checked" and b.improvements == []
    assert b.improvement == 0.0 and b.parent_total == pytest.approx(b.board_total)
    # removed a holder (1/5): regression
    with pytest.raises(VerifyError) as ei:
        verify_board({"r6.ec": R6}, parent, FAST, TABLE)
    assert ei.value.code is ErrorCode.BOARD_REGRESSION and "1/5" in ei.value.message
    # removed a breakthrough holder (1/6): regression
    with pytest.raises(VerifyError) as ei:
        verify_board({"b5.ec": B5}, parent, FAST, TABLE)
    assert ei.value.code is ErrorCode.BOARD_REGRESSION and "1/6" in ei.value.message
    # removed a non-holder: allowed
    r6_minus = "\n".join(l for l in R6.splitlines() if not l.startswith("P -633 10822")) + "\n"
    parent2 = {"b5.ec": B5, "r5big.ec": r6_minus}
    b = verify_board({"b5.ec": B5}, parent2, FAST, TABLE)
    assert b.parent_gate == "checked"


def test_improvements_are_contribution_deltas_against_the_parent():
    parent = {"g3.ec": G3}
    b = verify_board({"g3.ec": G3, "b5.ec": B5, "z2.ec": Z2R10}, parent, FAST, TABLE)
    imp = {i["cell"]: i for i in b.improvements}
    assert set(imp) == {"1/5", "2/10"}                       # 1/3: the parent's curve is the record itself, 0 before and after
    assert imp["1/5"]["file"] == "b5.ec" and imp["1/5"]["parent_contribution"] == 0.0
    assert imp["1/5"]["delta"] == pytest.approx(100 * (4.0 - 3.422446) / 4.0)
    assert imp["2/10"]["delta"] == pytest.approx(100 * (12.0 - 11.510206) / 12.0)
    assert b.parent_total == 0.0 and b.cells["1/3"].contribution == 0.0
    assert b.improvement == pytest.approx(b.board_total - b.parent_total)
    assert b.improvement == pytest.approx(imp["1/5"]["delta"] + imp["2/10"]["delta"])
    # a first curve in a cell neither the table nor the parent has: delta is the whole breakthrough term
    parent = {"b5.ec": B5}
    r6_rank4 = "\n".join(l for l in R6.splitlines() if not l.startswith(("P 2461 116018", "P -719/4"))) + "\n"
    b = verify_board({"b5.ec": B5, "b4.ec": r6_rank4}, parent, FAST, TABLE)
    assert [i["cell"] for i in b.improvements] == ["1/4"]
    assert b.improvements[0]["parent_contribution"] == 0.0
    assert b.improvements[0]["delta"] == pytest.approx(cells.contribution(4.288426, None)) == pytest.approx(b.improvement)


def test_parent_file_that_no_longer_verifies_rejects_the_run():
    # A promoted file verified when its frontier was accepted; with a frozen table and a pinned gp it
    # keeps verifying. One that fails now means the environment changed under the frontier (a
    # ledger rebuild, a verifier change). Its cell cannot be checked, so the gate fails closed
    # instead of listing the file and letting the candidate drop whatever it held.
    parent = {"b5.ec": B5, "junk.ec": "E 0 0 0 0 0\nP 1 1\n"}
    with pytest.raises(VerifyError) as ei:
        verify_board({"b5.ec": B5}, parent, FAST, TABLE)
    assert ei.value.code is ErrorCode.BOARD_REGRESSION
    assert "junk.ec" in ei.value.message and "CURVE_SINGULAR" in ei.value.message


def test_board_is_deterministic():
    a = verify_board({"b5.ec": B5, "r6.ec": R6}, None, FAST, TABLE).to_json()
    b = verify_board({"r6.ec": R6, "b5.ec": B5}, None, FAST, TABLE).to_json()
    for d in (a, b):
        for c in d["curves"]:
            c.pop("conductor")
    assert json.dumps(a, sort_keys=True) == json.dumps(b, sort_keys=True)


def test_rank_zero_curve_is_listed_but_contributes_nothing_and_is_never_held():
    nopts = "E 0 0 0 18372 -33031\n"           # base5's curve with no points: certified rank 0
    b = verify_board({"nopts.ec": nopts}, None, FAST, TABLE)
    assert b.board_total == 0.0
    assert b.cells["1/0"].contribution == 0.0
    js = b.to_json()
    (cur,) = js["curves"]
    assert cur["file"] == "nopts.ec" and cur["cell"] == "1/0" and cur["holds_cell"] is False
    assert js["cells"]["1/0"]["scored"] is False and js["cells"]["1/3"]["scored"] is True   # rank 0 is not a breakthrough
    by_file = {c["file"]: c for c in js["curves"]}
    assert by_file["nopts.ec"]["contribution"] == 0.0
    assert js["improvements"] == []
    b2 = verify_board({}, {"nopts.ec": nopts}, FAST, TABLE)
    assert b2.parent_gate == "checked" and b2.board_total == 0.0


def test_holder_that_gains_a_point_moves_up_without_regression():
    # A parent holder (rank6_blind minus two points: cell 1/4, a breakthrough) is resubmitted with
    # its two missing points (cell 1/6). Cell 1/4 is then empty on the candidate board, but the
    # curve that held it is still there at a higher rank; keeping the 4-point copy as well would
    # be DUPLICATE_CURVE, so this must not be a regression.
    r6_rank4 = "\n".join(l for l in R6.splitlines() if not l.startswith(("P 2461 116018", "P -719/4"))) + "\n"
    b = verify_board({"r6.ec": R6}, {"r4.ec": r6_rank4}, FAST, TABLE)
    assert b.parent_gate == "checked" and b.cells["1/6"].board_best_file == "r6.ec"
    # The reverse (the same curve with fewer certified points) is still a regression.
    with pytest.raises(VerifyError) as ei:
        verify_board({"r4.ec": r6_rank4}, {"r6.ec": R6}, FAST, TABLE)
    assert ei.value.code is ErrorCode.BOARD_REGRESSION and "1/6" in ei.value.message
    # Moving up is allowed only when the board total does not fall: the holder of breakthrough
    # cell 1/4 (worth about 192 under v6) resubmitted with one more point lands in 1/5 above that cell's
    # baseline (4.288 > 4.0), holds nothing there, and would vacate 1/4 for nothing.
    r6_rank5 = "\n".join(l for l in R6.splitlines() if not l.startswith("P -633 10822")) + "\n"
    with pytest.raises(VerifyError) as ei:
        verify_board({"r5.ec": r6_rank5}, {"r4.ec": r6_rank4}, FAST, TABLE)
    assert ei.value.code is ErrorCode.BOARD_REGRESSION and "1/4" in ei.value.message


def test_parent_file_matching_the_frozen_ledger_rejects_the_run():
    # The ledger is frozen (docs/records.md §4): a promoted curve can never become known later, so a parent
    # file that matches the ledger means the ledger was rebuilt under the frontier. That is a
    # maintainer's problem, surfaced as a rejection, not a cell silently unheld.
    with pytest.raises(VerifyError) as ei:
        verify_board({"b5.ec": B5}, {"g.ec": G3}, FAST_NOVELTY, TABLE)
    assert ei.value.code is ErrorCode.BOARD_REGRESSION and "KNOWN_CURVE" in ei.value.message


def test_holder_moving_up_is_credited_with_its_new_cells_gain_even_when_another_curve_is_best_there():
    # submitting.md §3: a holder may vacate its cell "as long as the board total does not fall". The
    # parent holds breakthrough 1/4 with rank6_blind minus two points. The candidate resubmits that
    # curve with all six points (cell 1/6) *and* adds the smaller known rank-6 record curve (novelty
    # is off under FAST) in the same cell. The moved curve is not the best in 1/6, but the board
    # total rose by the new breakthrough, so this is not a regression.
    import json, sys
    sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "tools"))
    from certify_known import submission_text
    raw = json.loads((pathlib.Path(__file__).resolve().parent.parent / "tools" / "known_raw.json").read_text())
    rec6 = submission_text(next(c for c in raw["curves"] if c["id"] == "icarm-56"))      # 1/6, H 2.531337
    r6_rank4 = "\n".join(l for l in R6.splitlines() if not l.startswith(("P 2461 116018", "P -719/4"))) + "\n"
    b = verify_board({"r6.ec": R6, "rec6.ec": rec6}, {"r4.ec": r6_rank4}, FAST, TABLE)
    assert b.cells["1/6"].board_best_file == "rec6.ec" and "1/4" not in b.cells     # vacated: not in the report
    assert b.board_total >= b.parent_total
    # If the smaller curve was already on the parent board, moving up gains nothing: regression.
    with pytest.raises(VerifyError) as ei:
        verify_board({"r6.ec": R6, "rec6.ec": rec6}, {"r4.ec": r6_rank4, "rec6.ec": rec6}, FAST, TABLE)
    assert ei.value.code is ErrorCode.BOARD_REGRESSION and "1/4" in ei.value.message


def test_two_holders_moving_into_one_cell_cannot_both_claim_its_gain():
    # Parent holds two breakthrough cells: 1/4 with rank6_blind minus two points (~192) and 1/2 with
    # the known rank-6 record curve minus four points (~195), total ~387. The candidate resubmits
    # both with all six points: both land in 1/6, worth ~195 once. Crediting that one gain to each
    # vacated cell would pass a board whose total fell to 195. The documented rule is global: the
    # total must not fall.
    import json, sys
    sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "tools"))
    from certify_known import submission_text
    raw = json.loads((pathlib.Path(__file__).resolve().parent.parent / "tools" / "known_raw.json").read_text())
    rec6 = submission_text(next(c for c in raw["curves"] if c["id"] == "icarm-56"))
    rec6_rank2 = "\n".join(rec6.splitlines()[:3]) + "\n"                 # E line + two points
    r6_rank4 = "\n".join(l for l in R6.splitlines() if not l.startswith(("P 2461 116018", "P -719/4"))) + "\n"
    parent = {"x4.ec": r6_rank4, "y2.ec": rec6_rank2}
    pb = verify_board(dict(parent), parent, FAST, TABLE)
    assert set(pb.cells) >= {"1/4", "1/2"} and pb.board_total > 380
    with pytest.raises(VerifyError) as ei:
        verify_board({"x6.ec": R6, "y6.ec": rec6}, parent, FAST, TABLE)
    assert ei.value.code is ErrorCode.BOARD_REGRESSION and "total" in ei.value.message
    # Moving only one of them up, keeping the other, is fine: 1/2 stays, 1/4 -> 1/6 gains 192 for 192.
    b = verify_board({"x6.ec": R6, "y2.ec": rec6_rank2}, parent, FAST, TABLE)
    assert b.board_total >= pb.board_total


def test_parent_second_pass_skips_the_conductor(monkeypatch):
    # The conductor is informational and wall-clock bounded; parent files that are only re-verified
    # for the gate must not spend the board budget on it. One call: the single candidate file.
    import ec_verify.verify as vmod
    calls = []
    monkeypatch.setattr(vmod.curve_mod, "conductor_best_effort", lambda *a, **k: calls.append(1) or None)
    r6_minus = "\n".join(l for l in R6.splitlines() if not l.startswith("P -633 10822")) + "\n"
    cfg = VerifyConfig(novelty_enabled=False)                 # conductor on, as in the official run
    verify_board({"b5.ec": B5}, {"b5.ec": B5, "r5big.ec": r6_minus}, cfg, TABLE)
    assert len(calls) == 1


@pytest.mark.parametrize("exc, code", [("GPTimeout", ErrorCode.RESOURCE_EXCEEDED), ("GPError", ErrorCode.GP_ERROR)])
def test_gp_failure_while_reverifying_a_parent_file_rejects_the_run(monkeypatch, exc, code):
    # A gp timeout or crash while re-verifying a parent holder is this run's environment, not a
    # parent file that "no longer verifies". Swallowing it into ignored_parent_files would unhold
    # the cell and let a candidate that dropped the holder pass with parent_gate "checked".
    import ec_verify.verify as vmod
    from ec_verify.parse import parse
    real = vmod.regulator.certify
    b5_ainvs = list(parse(B5).ainvs)

    def flaky(ainvs, free, gp_bin=None):
        if list(ainvs) == b5_ainvs:
            raise getattr(vmod.gpmod, exc)("gp exceeded 300s" if exc == "GPTimeout" else "GP_ERROR: gp crashed")
        return real(ainvs, free, gp_bin=gp_bin)

    monkeypatch.setattr(vmod.regulator, "certify", flaky)
    with pytest.raises(VerifyError) as ei:
        verify_board({"r6.ec": R6}, {"b5.ec": B5, "r6.ec": R6}, FAST, TABLE)
    assert ei.value.code is code and "b5.ec" in ei.value.message


def test_budget_is_checked_after_the_last_verification_too():
    # A board at 1199 s may start its last curve; the 1200 s cap must still reject the run when
    # that curve finishes late, not score it because no further tick was due.
    calls = []

    def check():
        calls.append(1)
        if len(calls) > 1:                      # the pre-verify tick passes; the one after must fire
            raise VerifyError(ErrorCode.RESOURCE_EXCEEDED, "board verification exceeded 1200 s")

    with pytest.raises(VerifyError) as ei:
        verify_board({"b5.ec": B5}, None, FAST, TABLE, budget_check=check)
    assert ei.value.code is ErrorCode.RESOURCE_EXCEEDED and len(calls) == 2


def test_board_run_does_not_relaunch_gp_for_the_version_per_curve(monkeypatch):
    # The harness checks the pin once; with the version in the config no curve pays another gp
    # launch for it (a 512-file board plus its parent would otherwise spend ~1000 launches of the
    # 1200 s budget on `gp --version-short`).
    import ec_verify.gp as gpm
    gp_bin = gpm.find_gp()
    launches = []
    for name in ("gp_version", "require_pinned"):
        monkeypatch.setattr(gpm, name, lambda b, _r=getattr(gpm, name): launches.append(1) or _r(b))
    cfg = VerifyConfig(conductor_enabled=False, novelty_enabled=False, gp_bin=gp_bin, pari_version=gpm.PINNED_VERSION)
    b = verify_board({"b5.ec": B5, "r6.ec": R6}, {"b5.ec": B5, "r6.ec": R6}, cfg, TABLE)
    assert launches == [] and b.pari_version == gpm.PINNED_VERSION
    assert all(r.pari_version == gpm.PINNED_VERSION for _, r in b.curves)


def test_verify_board_dirs_is_the_one_wiring_for_harness_and_cli(tmp_path, monkeypatch):
    # Strict loader for the board, lenient for the parent, a set parent that is not a directory is a
    # rejection (never an empty parent), the budget as a tick and a gp deadline. Both entry points
    # call this; neither hand-rolls it.
    import ec_verify.board as bmod
    board = tmp_path / "board"; board.mkdir(); (board / "b5.ec").write_text(B5, encoding="ascii")
    with pytest.raises(VerifyError) as ei:
        bmod.verify_board_dirs(board, tmp_path / "nope", FAST, TABLE)
    assert ei.value.code is ErrorCode.SUBMISSION_NOT_REGULAR_FILE and "nope" in ei.value.message
    b = bmod.verify_board_dirs(board, None, FAST, TABLE)
    assert b.parent_gate == "skipped" and [f for f, _ in b.curves] == ["b5.ec"]
    monkeypatch.setattr(bmod, "WALL_BUDGET_S", -1.0)
    with pytest.raises(VerifyError) as ei:
        bmod.verify_board_dirs(board, None, FAST, TABLE)
    assert ei.value.code is ErrorCode.RESOURCE_EXCEEDED


def test_environment_codes_are_defined_beside_the_error_codes():
    # The board layer must not keep its own list of which codes mean "this run's environment";
    # a code added to result.py's environment group must be picked up without touching board.py.
    import ec_verify.board as bmod
    from ec_verify import result
    assert bmod.ENVIRONMENT_CODES is result.ENVIRONMENT_CODES
    assert {ErrorCode.GP_UNAVAILABLE, ErrorCode.GP_VERSION_MISMATCH, ErrorCode.GP_ERROR, ErrorCode.RESOURCE_EXCEEDED} == set(result.ENVIRONMENT_CODES)


def test_empty_board_reports_the_configured_pari_version():
    cfg = VerifyConfig(novelty_enabled=False, conductor_enabled=False, pari_version="2.17.4")
    assert verify_board({}, None, cfg, TABLE).pari_version == "2.17.4"


def test_internal_error_in_a_candidate_is_gp_error_with_the_file_prefix(monkeypatch):
    # Garbage from gp (a KeyError/ValueError while reading its output) must reach the harness as the
    # documented GP_ERROR naming the file, not as a bare Python exception.
    from ec_verify import board as board_mod

    def boom(text, cfg):
        raise KeyError("AGREE")
    monkeypatch.setattr(board_mod, "verify_text", boom)
    with pytest.raises(VerifyError) as ei:
        verify_board({"a.ec": B5}, None, FAST, TABLE)
    assert ei.value.code is ErrorCode.GP_ERROR and ei.value.message.startswith("board/a.ec:") and "KeyError" in ei.value.message


def test_internal_error_in_a_parent_file_is_gp_error_naming_the_parent(monkeypatch):
    from ec_verify import board as board_mod

    def boom(text, cfg):
        raise ValueError("could not convert string to float: 'q'")
    monkeypatch.setattr(board_mod, "verify_text", boom)
    with pytest.raises(VerifyError) as ei:
        verify_board({}, {"p.ec": B5}, FAST, TABLE)
    assert ei.value.code is ErrorCode.GP_ERROR and ei.value.message.startswith("re-verifying parent board/p.ec:")


def test_identical_candidate_texts_are_verified_once(monkeypatch):
    from ec_verify import board as board_mod
    calls = []
    real = board_mod.verify_text

    def counting(text, cfg):
        calls.append(text)
        return real(text, cfg)
    monkeypatch.setattr(board_mod, "verify_text", counting)
    with pytest.raises(VerifyError) as ei:
        verify_board({"a.ec": B5, "a-copy.ec": B5}, None, FAST, TABLE)
    assert ei.value.code is ErrorCode.DUPLICATE_CURVE
    assert len(calls) == 1
