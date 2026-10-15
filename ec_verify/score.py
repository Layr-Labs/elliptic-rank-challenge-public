"""The Yukon scalar (docs/verifier-architecture.md §4): the board total. SCORING_VERSION changes whenever the formula does."""

from __future__ import annotations

SCORING_VERSION = 6


def yukon_score(board) -> float:
    """board is an ec_verify.board.BoardResult; the scalar is its total, unrounded."""
    return float(board.board_total)
