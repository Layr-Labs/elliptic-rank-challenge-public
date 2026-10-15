"""Error codes and the per-curve CurveResult record. Codes are API: agents parse them."""

from __future__ import annotations

import enum
import json
from dataclasses import dataclass


class ErrorCode(str, enum.Enum):
    # file / parse
    SUBMISSION_NOT_REGULAR_FILE = "SUBMISSION_NOT_REGULAR_FILE"
    PARSE_SYNTAX = "PARSE_SYNTAX"
    PARSE_EMPTY = "PARSE_EMPTY"
    CURVE_TOO_LARGE = "CURVE_TOO_LARGE"
    POINT_TOO_LARGE = "POINT_TOO_LARGE"
    TOO_MANY_POINTS = "TOO_MANY_POINTS"
    # curve
    CURVE_SINGULAR = "CURVE_SINGULAR"
    CURVE_NON_MINIMAL = "CURVE_NON_MINIMAL"
    KNOWN_CURVE = "KNOWN_CURVE"
    CURVE_TOO_MANY_BAD_PRIMES = "CURVE_TOO_MANY_BAD_PRIMES"
    # points
    POINT_OFF_CURVE = "POINT_OFF_CURVE"
    POINT_DUPLICATE = "POINT_DUPLICATE"
    POINTS_DEPENDENT = "POINTS_DEPENDENT"
    HEIGHT_PRECISION_DISAGREE = "HEIGHT_PRECISION_DISAGREE"
    # environment / resources
    GP_UNAVAILABLE = "GP_UNAVAILABLE"
    GP_VERSION_MISMATCH = "GP_VERSION_MISMATCH"
    GP_ERROR = "GP_ERROR"
    RESOURCE_EXCEEDED = "RESOURCE_EXCEEDED"
    # board
    DUPLICATE_CURVE = "DUPLICATE_CURVE"
    TOO_MANY_CURVES = "TOO_MANY_CURVES"
    BOARD_REGRESSION = "BOARD_REGRESSION"


# The "environment / resources" group above: codes that describe this run, not a file. The board
# layer re-raises these while re-verifying a parent file instead of treating the file as failed.
ENVIRONMENT_CODES = frozenset({ErrorCode.GP_UNAVAILABLE, ErrorCode.GP_VERSION_MISMATCH, ErrorCode.GP_ERROR,
                               ErrorCode.RESOURCE_EXCEEDED})


def internal_error_detail(e: BaseException) -> str:
    """The detail every entry point prints for an exception that is not a VerifyError (GP_ERROR,
    never a traceback): the exception type and the first 200 characters of its message."""
    return f"internal error: {type(e).__name__}: {str(e)[:200]}"


class VerifyError(Exception):
    def __init__(self, code: ErrorCode, message: str):
        super().__init__(f"{code.value}: {message}")
        self.code = code
        self.message = message

    def to_json(self) -> dict:
        return {"error": self.code.value, "message": self.message}


@dataclass(frozen=True)
class CurveResult:
    """Everything recorded per accepted curve."""

    ainvs: list                      # ints
    certified_rank: int
    torsion_order: int
    torsion: list                    # PARI elltors(E)[2], e.g. [6, 2]
    torsion_key: str                 # "6x2"
    cell: str                        # "6x2/4"
    naive_height: float              # H, 6 decimals
    naive_height_source: float       # 12*H (ICARM normalisation), 6 decimals
    naive_height_log10: float
    discriminant_digits: int
    fingerprint: list                # |a_p| for p <= 1000, always computed
    points_submitted: int
    points_torsion_dropped: int
    regulator: float
    regulator_lower_bound: float
    height_lambda_min: float | None
    height_precision_agreement: float
    claimed_rank: int | None
    claim_discrepancy: bool
    verified_claim: str
    conductor: str | None            # informational
    pari_version: str

    def to_json(self) -> dict:
        return {
            "ainvs": [str(a) for a in self.ainvs],
            "certified_rank": self.certified_rank,
            "torsion_order": self.torsion_order,
            "torsion": list(self.torsion),
            "torsion_key": self.torsion_key,
            "cell": self.cell,
            "naive_height": self.naive_height,
            "naive_height_source": self.naive_height_source,
            "naive_height_log10": self.naive_height_log10,
            "discriminant_digits": self.discriminant_digits,
            "fingerprint": list(self.fingerprint),
            "points_submitted": self.points_submitted,
            "points_torsion_dropped": self.points_torsion_dropped,
            "regulator": self.regulator,
            "regulator_lower_bound": self.regulator_lower_bound,
            "height_lambda_min": self.height_lambda_min,
            "height_precision_agreement": self.height_precision_agreement,
            "claimed_rank": self.claimed_rank,
            "claim_discrepancy": self.claim_discrepancy,
            "verified_claim": self.verified_claim,
            # wall-clock capped, machine dependent, never scored
            "conductor": self.conductor,
            "pari_version": self.pari_version,
        }

    def to_json_str(self) -> str:
        return json.dumps(self.to_json(), sort_keys=True, separators=(",", ":"))
