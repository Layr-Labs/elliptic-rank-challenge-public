"""The record-store pipeline tools (docs/records.md §2) on tiny fixture mirrors; gp-backed."""

import json
import pathlib
import sys

import pytest

from conftest import needs_gp
from util import base5
from ec_verify.parse import parse

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "tools"))
import certify_known  # noqa: E402

pytestmark = needs_gp
NOW = "2026-10-08T12:00:00+00:00"


def _entry(cid, ainvs, points, rank, torsion=None, updated="u1"):
    return {"id": cid, "source": cid.split("-")[0], "source_id": cid.split("-", 1)[1], "ainvs": ainvs, "points": points,
            "rank_claimed": rank, "torsion_source": torsion, "naive_height_source": None,
            "source_created_at": None, "source_updated_at": updated, "first_seen": NOW}


def _pts(text):
    return [[str(x) if x.denominator == 1 else f"{x.numerator}/{x.denominator}",
             str(y) if y.denominator == 1 else f"{y.numerator}/{y.denominator}"] for x, y in parse(text).points]


E5077 = _entry("lmfdb-5077a1", [0, 0, 1, -7, 6], [["-2", "3"], ["-1", "3"], ["0", "2"]], 3, [])
B5 = _entry("icarm-900001", list(parse(base5()).ainvs), _pts(base5()), 5, [])
NOPTS = _entry("lmfdb-11a1", [0, -1, 1, -10, -20], [], 0, [5])
DEP = _entry("icarm-900002", [0, 0, 1, -7, 6], [["-2", "3"], ["-1", "3"], ["0", "2"], ["1", "0"]], 4, [])   # (1,0) = -(P1+P2+... ) dependent


def test_submission_text_keeps_rational_points():
    txt = certify_known.submission_text(B5)
    assert txt.startswith("E 0 0 0 18372 -33031\n") and "P 764/25 93287/125\n" in txt


def test_certify_entry_statuses():
    gp = certify_known.gp_bin()
    v = certify_known.certify_entry(E5077, gp, NOW)
    assert v["status"] == "certified" and v["certified_rank"] == 3 and v["pari_version"] == "2.17.4"
    v = certify_known.certify_entry(NOPTS, gp, NOW)
    assert v["status"] == "unverified" and v["certified_rank"] is None
    v = certify_known.certify_entry(DEP, gp, NOW)
    assert v["status"] == "failed" and v["code"] == "POINTS_DEPENDENT" and v["certified_rank"] is None


def test_certify_all_is_incremental():
    raw = {"version": 2, "fetched_at": NOW, "curves": [dict(E5077), dict(NOPTS)]}
    n = certify_known.certify_all(raw, certify_known.gp_bin(), jobs=2, now=NOW)
    assert n == 2 and all("verification" in c for c in raw["curves"])
    assert certify_known.certify_all(raw, certify_known.gp_bin(), jobs=2, now=NOW) == 0
    raw["curves"][0]["source_updated_at"] = "u2"
    assert certify_known.certify_all(raw, certify_known.gp_bin(), jobs=2, now=NOW) == 1


import build_known  # noqa: E402
from ec_verify import known_curves


def _raw(*entries):
    return {"version": 2, "fetched_at": NOW, "curves": [dict(e) for e in entries]}


def test_build_ledger_fields_and_torsion_check(tmp_path):
    gp = certify_known.gp_bin()
    raw = _raw(E5077, NOPTS, B5)
    certify_known.certify_all(raw, gp, 2, NOW)
    ledger = build_known.build_ledger(raw, None, gp, jobs=2)
    known_curves.validate_ledger(ledger)
    by = {c["id"]: c for c in ledger["curves"]}
    assert by["lmfdb-11a1"]["torsion"] == [5] and by["lmfdb-11a1"]["torsion_key"] == "5"
    assert by["lmfdb-5077a1"]["verification"]["status"] == "certified" and by["lmfdb-5077a1"]["verification"]["certified_rank"] == 3
    assert "message" not in by["lmfdb-5077a1"]["verification"]
    assert by["icarm-900001"]["naive_height"] == pytest.approx(3.422446)
    assert [c["id"] for c in ledger["curves"]] == ["icarm-900001", "lmfdb-5077a1", "lmfdb-11a1"]   # -rank_claimed, id
    # a source torsion that disagrees with gp is fatal and names the id
    bad = _raw(dict(E5077, torsion_source=[2]))
    with pytest.raises(SystemExit) as ei:
        build_known.build_ledger(bad, None, gp, jobs=1)
    assert "lmfdb-5077a1" in str(ei.value)


def test_build_ledger_orders_a_dujella_copy_after_the_icarm_entry_of_the_same_curve():
    icarm = _entry("icarm-7", [0, 0, 1, -7, 6], [], 3, [])
    duj = _entry("dujella-z1-0a1b2c3d", [0, 0, 1, -7, 6], [], 3, None, updated=None)
    led = build_known.build_ledger(_raw(duj, icarm), None, certify_known.gp_bin(), 2)
    assert [c["id"] for c in led["curves"]] == ["icarm-7", "dujella-z1-0a1b2c3d"]
    assert [c["duplicate_of_earlier"] for c in led["curves"]] == [False, True]


def test_build_ledger_reuses_previous_fingerprints(tmp_path, monkeypatch):
    gp = certify_known.gp_bin()
    raw = _raw(E5077, B5)
    certify_known.certify_all(raw, gp, 2, NOW)
    first = build_known.build_ledger(raw, None, gp, jobs=2)
    calls = []
    import ec_verify.fingerprint as fpm
    real = fpm.fingerprint
    monkeypatch.setattr(build_known.fp_mod, "fingerprint", lambda a, gp_bin=None: calls.append(a) or real(a, gp_bin=gp_bin))
    second = build_known.build_ledger(raw, first, gp, jobs=2)
    assert calls == []                                        # both fingerprints reused from `first`
    by = {c["id"]: c for c in second["curves"]}
    assert "novelty_exempt" not in by["icarm-900001"]
    assert known_curves.match(by["icarm-900001"]["fingerprint"], second)["id"] == "icarm-900001"


import build_records  # noqa: E402
from ec_verify import records as records_mod


def _ledger(*curves):
    return {"version": 2, "frozen_at": NOW, "primes_upto": 1000, "curves": list(curves)}


def _lc(cid, tkey, torsion, rank, H, fp, status="certified"):
    return {"id": cid, "source": cid.split("-")[0], "source_id": cid.split("-", 1)[1], "ainvs": ["0", "0", "0", "1", str(H)],
            "fingerprint": fp, "torsion": torsion, "torsion_key": tkey, "naive_height": H, "rank_claimed": rank,
            "verification": {"status": status, "certified_rank": rank if status == "certified" else None, "code": None,
                             "verified_at": NOW, "pari_version": "2.17.4", "verifier_version": "x"},
            "duplicate_of_earlier": False, "source_created_at": None, "first_seen": NOW}


FPA, FPB, FPC = [1] * 168, [2] * 168, [3] * 168


def test_build_table_records_counts_and_validation():
    led = _ledger(_lc("icarm-1", "2", [2], 9, 13.0, FPA), _lc("icarm-2", "2", [2], 9, 12.5, FPB),
                  _lc("icarm-3", "2", [2], 9, 12.5, FPB),                      # isogenous to icarm-2: same class
                  _lc("icarm-4", "1", [], 6, 4.0, FPC, status="failed"))      # never defines a record
    t = build_records.build_table(led, None, NOW)
    records_mod.validate(t)
    assert set(t["cells"]) == {f"2/{r}" for r in range(1, 10)} and t["version"] == 1   # 2/1..2/8 filled from 2/9
    c = t["cells"]["2/9"]
    assert c["record"]["id"] == "icarm-2" and c["record"]["certified_rank"] == 9 and c["exact_count"] == 3
    assert c["known_count"] == 3 and c["class_count"] == 2
    assert "baseline" not in c and "h_cap" not in t


def test_every_rank_up_to_the_top_is_a_cell_whose_record_is_the_smallest_curve_reaching_it():
    # CLAUDE.md: a cell is a breakthrough when *no known curve of that torsion reaches the rank*.
    # A certified rank-11 Z/2 curve reaches rank 10 (list ten of its points), so 2/10 must exist,
    # with the smallest curve of rank >= 10 as its record; 2/12 and above stay breakthroughs. The
    # same holds below the lowest exact cell: 2/1..2/8 take the rank-9 curve.
    led = _ledger(_lc("icarm-1", "2", [2], 9, 13.0, FPA), _lc("icarm-2", "2", [2], 11, 20.0, FPB),
                  _lc("icarm-3", "2", [2], 11, 18.0, FPC))
    t = build_records.build_table(led, None, NOW)
    records_mod.validate(t)
    assert set(t["cells"]) == {f"2/{r}" for r in range(1, 12)}
    low = t["cells"]["2/3"]
    assert low["record"]["id"] == "icarm-1" and low["record"]["certified_rank"] == 9
    assert low["exact_count"] == 0 and low["known_count"] == 3
    gap = t["cells"]["2/10"]
    assert gap["rank"] == 10 and gap["record"]["id"] == "icarm-3" and gap["record"]["naive_height"] == 18.0
    assert gap["known_count"] == 2 and gap["class_count"] == 2 and gap["exact_count"] == 0
    assert gap["record"]["certified_rank"] == 11 and t["cells"]["2/11"]["exact_count"] == 2
    # A curve certified at exactly rank 10 but larger (H 25) than the rank-11 curve (H 18) does
    # not become the record: the rank-11 curve still reaches rank 10 and is smaller. (The v5
    # tool took the exact-rank minimum here; two shipped cells had a baseline above a known curve.)
    led2 = _ledger(*led["curves"], _lc("icarm-4", "2", [2], 10, 25.0, FPA))
    t2 = build_records.build_table(led2, t, NOW)
    c210 = t2["cells"]["2/10"]
    assert c210["record"]["id"] == "icarm-3" and c210["exact_count"] == 1 and c210["known_count"] == 3
    # ... and a smaller exact-rank curve does
    led3 = _ledger(*led["curves"], _lc("icarm-5", "2", [2], 10, 15.0, FPB))
    assert build_records.build_table(led3, t, NOW)["cells"]["2/10"]["record"]["id"] == "icarm-5"


def test_build_table_breaks_a_height_tie_by_source_then_numeric_id():
    # The same curve mirrored from ICARM and from a Dujella page has one height: the record keeps
    # the ICARM id (so a refresh that adds Dujella copies does not reshuffle the table), and among
    # ICARM copies the numerically smallest id wins (icarm-9 before icarm-10).
    led = _ledger(_lc("dujella-z2-0a1b2c3d", "2", [2], 9, 12.5, FPA), _lc("icarm-10", "2", [2], 9, 12.5, FPA),
                  _lc("icarm-9", "2", [2], 9, 12.5, FPA), _lc("lmfdb-37a1", "2", [2], 9, 12.5, FPA))
    assert build_records.build_table(led, None, NOW)["cells"]["2/9"]["record"]["id"] == "icarm-9"
    led = _ledger(_lc("dujella-z2-0a1b2c3d", "2", [2], 9, 12.5, FPA), _lc("lmfdb-37a1", "2", [2], 9, 12.5, FPA))
    assert build_records.build_table(led, None, NOW)["cells"]["2/9"]["record"]["id"] == "dujella-z2-0a1b2c3d"


def test_build_table_never_builds_a_rank_zero_cell():
    led = _ledger(_lc("icarm-1", "1", [], 0, 3.0, FPA), _lc("icarm-2", "2", [2], 9, 12.5, FPB))
    t = build_records.build_table(led, None, NOW)
    assert "1/0" not in t["cells"] and not any(c["torsion_key"] == "1" for c in t["cells"].values())
    assert {k for k in t["cells"] if k.startswith("2/")} == {f"2/{r}" for r in range(1, 10)}


# --- serialisation: one record per line so refresh diffs stay reviewable -------------------

def test_dump_per_line_writes_one_record_per_line_and_round_trips(tmp_path):
    import fetch_known as fk
    doc = {"version": 2, "fetched_at": "t", "curves": [{"id": "b", "ainvs": [1, 2]}, {"id": "a", "points": [["1/2", "3"]]}]}
    text = fk.dump_per_line(doc, "curves")
    assert json.loads(text) == doc
    lines = text.splitlines()
    assert text.endswith("\n") and lines[0] == "{" and lines[-1] == "}"
    # one line per curve, compact, keys sorted; input order preserved
    curve_lines = [l for l in lines if l.startswith("{\"")]
    assert curve_lines == ['{"ainvs":[1,2],"id":"b"},', '{"id":"a","points":[["1/2","3"]]}']
    assert fk.dump_per_line(doc, "curves") == text                     # deterministic


def test_shipped_data_files_are_one_record_per_line():
    for rel in ("ec_verify/known_curves.json", "tools/known_raw.json"):
        text = (ROOT / rel).read_text()
        n = len(json.loads(text)["curves"])
        assert len(text.splitlines()) < n + 20, f"{rel}: {len(text.splitlines())} lines for {n} records"


def test_certify_all_reads_the_gp_version_once(monkeypatch):
    # One gp launch for the version string per run, not one per curve (5586 on a --force pass).
    class Counting:
        def __init__(self, mod):
            self._mod, self.calls = mod, 0

        def gp_version(self, gp):
            self.calls += 1
            return self._mod.gp_version(gp)

        def __getattr__(self, name):
            return getattr(self._mod, name)

    gp = certify_known.gp_bin()
    proxy = Counting(certify_known.gpmod)
    monkeypatch.setattr(certify_known, "gpmod", proxy)
    # ... and the per-curve verifier must not launch gp for the version either: count at the source.
    import ec_verify.gp as gpm
    launches = []
    for name in ("gp_version", "require_pinned"):
        monkeypatch.setattr(gpm, name, lambda b, _n=name, _r=getattr(gpm, name): launches.append(_n) or _r(b))
    raw = {"curves": [dict(E5077), dict(NOPTS), dict(DEP)]}
    assert certify_known.certify_all(raw, gp, jobs=1, now=NOW) == 3
    assert proxy.calls == 1
    assert len(launches) == 1, launches
    assert {c["verification"]["pari_version"] for c in raw["curves"]} == {certify_known.gpmod.PINNED_VERSION}


def test_build_records_main_takes_the_version_from_the_shipped_table_not_the_output_path(tmp_path, monkeypatch):
    # A preview build into another path must still number itself after the shipped table, or
    # records_version in score.json can go backwards when the preview is copied in.
    import json
    shipped = json.loads(build_records.OUT.read_text())
    out = tmp_path / "preview.json"
    monkeypatch.setattr("sys.argv", ["build_records.py", "--out", str(out)])
    assert build_records.main() == 0
    assert json.loads(out.read_text())["version"] == shipped["version"] + 1


def test_certify_entry_does_not_record_an_environment_failure_as_a_verdict(monkeypatch):
    # A gp timeout or crash while certifying a known curve says nothing about the curve; recording it
    # as `failed` would drop the curve from the table (and never retry). The run must stop instead.
    from ec_verify import result as rmod
    monkeypatch.setattr(certify_known, "verify_text",
                        lambda *a, **k: (_ for _ in ()).throw(rmod.VerifyError(rmod.ErrorCode.RESOURCE_EXCEEDED, "gp exceeded 300s")))
    with pytest.raises(rmod.VerifyError) as ei:
        certify_known.certify_entry(dict(E5077), certify_known.gp_bin(), NOW, pari_version="2.17.4")
    assert ei.value.code is rmod.ErrorCode.RESOURCE_EXCEEDED
