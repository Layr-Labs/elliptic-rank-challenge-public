"""Per-curve pipeline: text in, CurveResult out, VerifyError on any gate failure (docs/verifier-architecture.md §2).
Board-level logic (many files, duplicates, regression, total) is ec_verify/board.py."""

from __future__ import annotations

from dataclasses import dataclass

from . import cells, curve as curve_mod, fingerprint as fp_mod, gp as gpmod, known_curves
from . import points as points_mod, regulator, torsion as torsion_mod
from .parse import parse
from .result import CurveResult, ErrorCode, VerifyError


@dataclass(frozen=True)
class VerifyConfig:
    gp_bin: str | None = None
    require_pinned_gp: bool = True
    conductor_enabled: bool = True
    novelty_enabled: bool = True
    pari_version: str | None = None  # set by a caller that has already checked the pin: no per-curve gp launch for it
    deadline: float | None = None    # time.monotonic() instant after which no gp stage may run (the board budget)


def gp_guard(fn, *args, **kwargs):
    """Run a gp-backed step, mapping GPError subclasses to the documented environment codes."""
    try:
        return fn(*args, **kwargs)
    except gpmod.GPTimeout as e:
        raise VerifyError(ErrorCode.RESOURCE_EXCEEDED, str(e))
    except gpmod.GPVersionMismatch as e:
        raise VerifyError(ErrorCode.GP_VERSION_MISMATCH, str(e).split(": ", 1)[-1])
    except gpmod.GPError as e:
        msg = str(e)
        code = ErrorCode.GP_UNAVAILABLE if msg.startswith("GP_UNAVAILABLE") else ErrorCode.GP_ERROR
        raise VerifyError(code, msg.split(": ", 1)[-1])


def verify_text(text: str, config: VerifyConfig | None = None) -> CurveResult:
    cfg = config or VerifyConfig()
    with gpmod.deadline(cfg.deadline):
        return _verify_text(text, cfg)


def _verify_text(text: str, cfg: VerifyConfig) -> CurveResult:
    sub = parse(text)
    ainvs = sub.ainvs

    disc = curve_mod.require_nonsingular(ainvs)
    points_mod.check_points(ainvs, sub.points, sub.point_lines)

    gp_bin = cfg.gp_bin or gp_guard(gpmod.find_gp)
    if cfg.pari_version is not None:
        pari_version = cfg.pari_version
    elif cfg.require_pinned_gp:
        pari_version = gp_guard(gpmod.require_pinned, gp_bin)
    else:
        pari_version = gp_guard(gpmod.gp_version, gp_bin)

    gp_guard(curve_mod.require_minimal_at_small_primes, ainvs, disc, gp_bin=gp_bin)
    # The fingerprint is always computed: the board layer needs it for DUPLICATE_CURVE.
    fp = gp_guard(fp_mod.fingerprint, ainvs, gp_bin=gp_bin)
    # The good-prime floor is what makes fingerprint matching decidable, for the novelty gate and
    # for DUPLICATE_CURVE alike, so it holds whether or not the known-curve lookup runs.
    good = fp_mod.good_primes(fp)
    if good < fp_mod.CANDIDATE_MIN_GOOD:
        raise VerifyError(ErrorCode.CURVE_TOO_MANY_BAD_PRIMES,
                          f"good reduction at only {good} of the {len(fp_mod.PRIMES)} primes below "
                          f"{fp_mod.PRIME_BOUND}; at least {fp_mod.CANDIDATE_MIN_GOOD} are required for "
                          "fingerprint matching to be decidable (docs/submitting.md §4)")
    if cfg.novelty_enabled:
        hit = known_curves.match(fp)
        if hit is not None:
            raise VerifyError(ErrorCode.KNOWN_CURVE,
                              f"curve is {hit['id']} (claimed rank >= {hit['rank_claimed']}) or isogenous/twist-equivalent "
                              "to it; known curves are not eligible (docs/submitting.md §4)")
    tors, structure, free, dropped = gp_guard(torsion_mod.strip_torsion, ainvs, list(sub.points), gp_bin=gp_bin)
    cert = gp_guard(regulator.certify, ainvs, free, gp_bin=gp_bin)
    rank = cert["rank"]
    try:
        tkey = cells.torsion_key(structure)
    except ValueError as e:                    # a STRUCT line gp should never produce: environment, not the curve
        raise VerifyError(ErrorCode.GP_ERROR, f"torsion structure from gp: {e}")
    cell = cells.cell_key(tkey, rank)
    H = curve_mod.naive_height(ainvs)
    N = curve_mod.conductor_best_effort(ainvs, disc, gp_bin=gp_bin) if cfg.conductor_enabled else None

    claim = f"certified rank >= {rank} from {len(free)} independent points of infinite order"
    if dropped:
        claim += f" ({len(dropped)} torsion point(s) not counted)"
    claim += f"; torsion {tkey}; cell {cell}; H = {H:.6f}"
    discrepancy = sub.claimed_rank is not None and sub.claimed_rank != rank
    if discrepancy:
        claim += f"; #CLAIM rank={sub.claimed_rank} not established"

    return CurveResult(
        ainvs=list(ainvs), certified_rank=rank, torsion_order=tors, torsion=list(structure), torsion_key=tkey,
        cell=cell, naive_height=H, naive_height_source=round(12.0 * H, 6),
        naive_height_log10=curve_mod.naive_height_log10(ainvs), discriminant_digits=curve_mod.digits(disc),
        fingerprint=fp, points_submitted=len(sub.points), points_torsion_dropped=len(dropped),
        regulator=cert["regulator"], regulator_lower_bound=cert["regulator_lower_bound"],
        height_lambda_min=cert["lambda_min"], height_precision_agreement=cert["agreement"],
        claimed_rank=sub.claimed_rank, claim_discrepancy=discrepancy, verified_claim=claim,
        conductor=N, pari_version=pari_version,
    )
