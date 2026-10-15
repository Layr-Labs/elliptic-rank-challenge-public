"""Elliptic Rank Challenge verifier.

Certifies each submitted curve (exact on-curve checks, torsion stripping, a
certified positive-definite canonical height pairing) and scores the board of
certified curves against the per-torsion, per-rank size records. Nothing in a
submission is trusted; everything is recomputed.
"""

from .result import CurveResult, ErrorCode, VerifyError

__version__ = "1.0.0.dev0"


def verify_text(text, config=None):
    """Lazy import so tools that only need gp load fast."""
    from .verify import verify_text as _vt
    return _vt(text, config)


def __getattr__(name):
    """`VerifyConfig` is exported lazily for the same reason."""
    if name == "VerifyConfig":
        from .verify import VerifyConfig
        return VerifyConfig
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


__all__ = ["CurveResult", "ErrorCode", "VerifyConfig", "VerifyError", "verify_text", "__version__"]
