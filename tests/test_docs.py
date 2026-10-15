"""Numbers the docs quote must be the code's."""

import pathlib
import re

from ec_verify import cells
from ec_verify.score import SCORING_VERSION

ROOT = pathlib.Path(__file__).resolve().parent.parent


def test_docs_quote_the_live_constants():
    for name in ("README.md", "AGENTS.md", "docs/submitting.md", "docs/verifier-architecture.md"):
        text = (ROOT / name).read_text()
        assert "f_mn" not in text and "calibration.json" not in text and "--fast" not in text, name
        assert "worth 100." not in text and "worth 100 " not in text, name     # v6: a breakthrough is worth more than 100
        assert f"{cells.MAX_BOARD_FILES}" in text, name
    arch = (ROOT / "docs/verifier-architecture.md").read_text()
    assert f"H_REF = {int(cells.H_REF)}" in arch and f"scoring version {SCORING_VERSION}" in arch


def test_benchmark_description_states_the_breakthrough_rule():
    import json
    desc = json.loads((ROOT / "benchmark.json").read_text())["description"]
    assert "more than 100" in desc and "or 100 for" not in desc


def test_agents_and_claude_are_identical():
    assert (ROOT / "AGENTS.md").read_bytes() == (ROOT / "CLAUDE.md").read_bytes()
