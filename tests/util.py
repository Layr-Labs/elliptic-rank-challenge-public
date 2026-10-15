"""Shared helpers: fixture texts and a fast VerifyConfig for gate tests.

Negative cases mutate a valid fixture so each isolates one property. No test
reads the live submission/ (participant-owned; see test_suite_hygiene)."""

from __future__ import annotations

import pathlib

from ec_verify.verify import VerifyConfig

ROOT = pathlib.Path(__file__).resolve().parent.parent
FIXTURES = ROOT / "tests" / "fixtures"
CORPUS = ROOT / "tests" / "corpus"

# Conductor off (wall-clock bounded extra). Gates and certification are identical to the official run.
FAST = VerifyConfig(conductor_enabled=False, novelty_enabled=False)
FAST_NOVELTY = VerifyConfig(conductor_enabled=False)


def fixture(name: str) -> str:
    return (FIXTURES / name).read_text(encoding="ascii")


def corpus(name: str) -> str:
    return (CORPUS / name).read_text(encoding="ascii")


def gens_5077a() -> str:
    return fixture("5077a_gens.ec")


def dependent_5077a() -> str:
    return fixture("5077a_dependent.ec")


def base5_dependent() -> str:
    return fixture("base5_dependent.ec")


def base5() -> str:
    """The novel rank-5 baseline curve (not in the frozen known set)."""
    return fixture("base5.ec")


def count_points(text: str) -> int:
    return sum(1 for l in text.splitlines() if l.startswith("P "))
