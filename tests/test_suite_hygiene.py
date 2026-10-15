"""No test may read the live submission/ — it is participant-owned and changes with every promotion."""

import pathlib
import re

HERE = pathlib.Path(__file__).resolve().parent


def test_no_test_reads_live_submission():
    for p in HERE.glob("test_*.py"):
        if p.name == "test_suite_hygiene.py":
            continue
        src = p.read_text()
        assert "submission/board" not in src and 'ROOT / "submission"' not in src, p.name


def test_fixtures_are_ascii_lf():
    for p in list((HERE / "fixtures").glob("*.ec")) + list((HERE / "corpus").glob("*.ec")):
        raw = p.read_bytes()
        assert raw.isascii() and b"\r" not in raw, p.name
