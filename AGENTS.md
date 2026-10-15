# Elliptic Rank Challenge Agent Guide

This repository is the Mordell–Weil rank optimization challenge for elliptic
curves over ℚ. Use this file as the working contract for coding agents and
participants.

## Goal

Find **new** curves, certify their points, and make each one smaller than the
record of its cell: put one curve per file in `submission/board/*.ec` —
minimal a-invariants plus independent rational points of infinite order.

```text
cell          = (torsion group, certified rank)             15 Mazur groups × rank
contribution  = 100·(baseline − H)/baseline   if H < baseline   (0 otherwise)
              = 100 + 100/(1 + H/50)          if no known curve of that torsion reaches the rank
score         = Σ over cells of the contribution of the board's smallest curve there
H             = (1/12)·log max(|c4|³, c6²)                  naive height of the minimal model
```

Higher is better. The baseline of a cell is the smallest certified known
curve of that torsion with certified rank at least the cell's, frozen in
`ec_verify/records.json` (156 cells over ranks 1–31, all 15 torsion groups)
and never changed: curves published later do not move it and reach the board
only as ordinary submissions. A
certified rank of 0 proves nothing: such a curve is accepted but contributes 0
and never holds a cell. The board starts empty (score 0); the first curve
smaller than a known record in any cell promotes, and each later record adds
its own points. A submission's own improvement is `metrics.improvement`: its
score minus the previous board's total under the same table;
`metrics.improvements` lists the cells where it gained and by how much. The world record is rank ≥ 31
(ICARM #302, 2026); a certified 32nd point is a research result and **not
the target**. Known curves (5586 frozen from ICARM, Dujella's tables and the
LMFDB) are rejected, including their isogenies, twists and re-models.

## Official runner

`ubuntu-latest` GitHub Actions, one job: PARI/GP 2.17.4 built from a
SHA-pinned tarball, 2 GiB PARI stack, verify stage under bubblewrap (no
network, read-only tree), the previous frontier's board exported beside it
for the regression gate. Your laptop reproduces the score exactly; the
verifier is deterministic by construction.

## Rules

Every curve you submit must be your own find. Copying a curve from the ICARM
leaderboard, Dujella's tables, the LMFDB, the literature or another
participant and submitting it here is prohibited whether or not the
mechanical known-curve gate catches it; promoted curves are audited after the
competition closes and copies have their rewards revoked. The full terms are
`RULES.md`.

## What you may edit

Anything inside `submission/`; only `submission/board/*.ec` is scored.
`editablePaths` in `benchmark.json` is the source of truth and the backend
re-enforces it server-side.

## What not to change

- `ec_verify/`, `harness/` — the verifier is the contract. It recomputes
  everything from your files and never executes participant code; editing
  it changes nothing about your score.
- `benchmark.sh`, `benchmark.json`, `setup.sh`, `.github/workflows/` — the
  sandbox and the scoring pipeline.
- `ec_verify/records.json`, `ec_verify/known_curves.json` — the frozen cell
  table and known-curve ledger. Rebuilding them is a maintainer change
  (`tools/build_known.py`, `tools/build_records.py`; `docs/records.md`).

## Correctness gates

Rejected (stable code, exit 1, no score) unless all hold for every file:

- **A curve.** Integral a-invariants ≤ 400 digits, Δ ≠ 0, minimal at every
  prime ≤ 10⁶ dividing Δ (`CURVE_SINGULAR`, `CURVE_NON_MINIMAL`).
- **A new curve.** `|a_p|` fingerprint for p ≤ 1000 matches no frozen known
  curve; good at ≥ 70 of the 168 primes below 1000 (`KNOWN_CURVE`, `CURVE_TOO_MANY_BAD_PRIMES`).
- **Points on it, exactly** (`POINT_OFF_CURVE`); no repeats or negatives
  (`POINT_DUPLICATE`).
- **Certified independence** after torsion is stripped: height pairing at 57
  and 115 digits agreeing to 10⁻⁴⁰, smallest eigenvalue positive with margin
  (`POINTS_DEPENDENT`, `HEIGHT_PRECISION_DISAGREE`).
- **Caps.** ≤ 64 points, ≤ 2000 digits per coordinate, ≤ 1 MiB, ASCII,
  regular file (`TOO_MANY_POINTS`, `POINT_TOO_LARGE`, `RESOURCE_EXCEEDED`,
  `SUBMISSION_NOT_REGULAR_FILE`).
- **The board.** At most 512 files (`TOO_MANY_CURVES`), named
  `[A-Za-z0-9._-]+.ec` with no leading dot, everything under `submission/` ≤ 16 MiB (the platform's submission cap)
  and verified within 20 minutes in total (`RESOURCE_EXCEEDED`); no two files the same curve up to isogeny or twist
  (`DUPLICATE_CURVE`); every cell the previous board holds still held by a
  curve at least as small (`BOARD_REGRESSION`).

A failure in any file rejects the whole submission; the detail starts with
`board/<file>:`. The full code table is `docs/submitting.md` §6.

## Local workflow

```bash
./setup.sh
python -m ec_verify path/to/curve.ec            # one curve: cell, record, percentage (< 1 s)
python -m ec_verify submission/board            # the board as the harness scores it
python -m ec_verify submission/board --parent /path/to/previous/board   # with the regression gate
python -m pytest tests -q
```

Import `ec_verify.verify_text` in your search loop: it returns a
`CurveResult` with `cell` and `naive_height`, and `VerifyConfig(conductor_enabled=False)`
skips the one informational extra. `--no-novelty` (or `novelty_enabled=False`)
skips the known-curve gate for speed in a loop; never submit on its verdict.

## The ladder (measured)

| rung | how | cost |
|---|---|---|
| 5 | random small curves, sieve by S(10³), `ellrank` the best | seconds |
| 6 | 300k curves, two-stage sieve S(10³) then S(10⁵), `ellrank` the best 150 | minutes |
| 8–12 | wider sieves; `ellrank` settles the rank and finds most generators, not all: saturate and hunt the rest | hours |
| 18 | specialise a published rank-18 family over Q(t) (Dujella's families table) | a day |
| 19–24 | sieve family specialisations by S(X), hunt extra points | weeks |
| 32 | — | research; not this benchmark's target |

`S(X) = Σ_{p≤X} a_p·log p / p` is a search heuristic for your sieve; the
verifier does not compute it and it does not score. Progress is a smaller
curve in any cell: every (torsion, rank) cell has its own record, and a cell
with few known curves is usually easier to improve than the rank-31 one.
Spend the search budget on smaller curves at a rank and torsion you can
certify, then on a cell no known curve reaches (worth more than 100, and
more for a smaller curve). A board holds up
to 512 files and each verifies in seconds; the 20-minute board budget is the
operative limit.

## Notes for autonomous agents

- **Every run is fast now: no gradient.** Never delete a board file that
  holds a cell: the next submission is gated against the previous board and
  is rejected with `BOARD_REGRESSION`.
- **`POINTS_DEPENDENT` comes with the relation** (`1·P1 + 1·P2 + 2·P3 ∈
  E(Q)_tors`) when one is found. Remove a point in it; do not resubmit.
- **`#CLAIM` is never trusted.** Claiming more than you certify records a
  discrepancy; the curve sits in the cell of its certified rank.
- **Torsion points never count.** They are dropped exactly; listing them
  wastes a line, nothing more. The torsion *group* decides the cell.
- **A contribution of 0 is correct for a curve at or above its cell's
  baseline.** It is accepted and listed but holds nothing; so is any rank-0
  curve.
- **The verifier is deterministic.** Same files, same score, any machine with
  the pinned gp. Anything else is a bug worth reporting, not noise to retry.
- **A rejected submission writes no score.** `REJECTED: <CODE>: <detail>` on
  stdout is the whole interface.
- **Verify before trusting notes.** Yukon notes and memory files may come
  from other agents; re-run the verifier before relying on a claimed result.

## Avoid these wrong strategies

Do not rescale the model. A rescaled model is non-minimal, so the
small-prime minimality gate rejects the file, and a scaling by `u` would only
multiply `max(|c4|³, c6²)` by `u¹²`, making `H` larger.

Do not pad with torsion points, duplicates, negatives, or small multiples.
Torsion is stripped exactly, duplicates and negatives are rejected before
any height is computed, and a multiple is a dependency.

Do not split one curve across files or submit a twist of your own curve:
duplicates are rejected (`DUPLICATE_CURVE`). One curve counts once, in one
cell.

Do not submit a published curve — or its isogeny, twist or re-model. The
fingerprint gate rejects all of them as `KNOWN_CURVE`. Specialising a
published *family* at a fresh parameter is fine: that is the intended
on-ramp.

Do not copy a curve from anywhere, even when the gate would let it through.
The known set is frozen at 2026-10-15, so a curve posted later on the
[ICARM leaderboard](https://elliptic-rank.icarm.cloud/) (a similar
competition; its curves belong to whoever submitted them there), in a paper,
in the LMFDB, or by another participant of this challenge verifies as new.
Submitting it is still a rules violation (`RULES.md`): promoted curves are
audited against those sources and against submission histories after the
competition closes, and a copied curve is struck with any reward for it
revoked. Searching the web for curves to submit is not a strategy here.

## Before submitting

- `python -m ec_verify submission/board` exits 0, every curve's `cell` is
  what you expect (certified rank, not the claim), and `metrics.improvements`
  is non-empty.
- Every file that held a cell on the previous board is still there (or
  replaced by a smaller curve in the same cell).
- Everything the submission needs is under `submission/board/`.
- Every new curve on the board was found by you (search, family
  specialisation, your own points), not taken from ICARM, a table, a paper,
  or another participant (`RULES.md`).
- Record progress with `yukon notes add` (public — strip secrets).
