# The Elliptic Rank Challenge

**Find a new curve. Prove its points. Beat the size record of its cell.**

> **Goal.** Submit elliptic curves over ℚ **nobody has published** together
> with independent rational points; the number of points the verifier can
> certify independent is a certified lower bound on the **Mordell–Weil rank**.
> Each curve lands in the cell of its torsion group and certified rank, and
> scores by the percentage by which it is **smaller** than the smallest
> certified known curve of that cell. A certified rank that no known curve of
> that torsion reaches is worth more than 100. Higher is better. The rank record (≥ 31)
> is context, not the target: the progress this challenge is built to measure
> is smaller curves at the ranks where new curves can actually be found and
> certified.
>
> **Rules.** Every curve you submit must be your own find. Copying a curve
> from anywhere — the [ICARM leaderboard](https://elliptic-rank.icarm.cloud/),
> Dujella's tables, the LMFDB, a paper, or another participant — and
> submitting it here is not allowed, whether or not the verifier's frozen
> known-curve gate catches it. Promoted curves are audited after the
> competition closes and copies have their rewards revoked. Full terms:
> [`RULES.md`](RULES.md).

---

## Why this matters

By Mordell's theorem `E(ℚ) ≅ E(ℚ)_tors × ℤ^r`; the integer `r` is the rank.
Whether ranks are bounded is open. Every record is a **lower bound** proven
by exhibiting `r` independent points — so the record *is* a witness, and a
machine can check it: the points lie on the curve exactly, none has finite
order, and the canonical height pairing of the rest is positive definite.

The record is **rank ≥ 31**: curve #302 on the
[ICARM Elliptic Curve Rank Leaderboard](https://elliptic-rank.icarm.cloud/curve/302),
posted 23 August 2026 and credited to Claude, Levent Alpöge and Ava Howell
(30 in the same month; 29 by Elkies–Klagsbrun 2024; 28 by Elkies 2006 — see
Andrej Dujella's [history of rank records](https://web.math.pmf.unizg.hr/~duje/tors/rankhist.html)).
This verifier certifies the 31 published points independent in under a
second, with no GRH, BSD or descent assumption; the one numerical input is
PARI's canonical height, cross-checked at two precisions
([`docs/verifier-architecture.md`](docs/verifier-architecture.md) §3 states
exactly what the certificate relies on). The largest rank known *exactly*,
unconditionally, is 22 (Samuel Li, September 2026, after Li's 21 and
Elkies–Klagsbrun's 20 of 2020).

**A 32nd point is not the target.** Each recent record came from a new
family found by a research programme, and one more such curve would move the
field very little: the open question is whether ranks are bounded, and a
single record does not answer it. Where progress is plausible, and where the
community already competes, is at the ranks below and in the torsion groups
other than the trivial one: for a given torsion group and rank, curves with
smaller coefficients and better-chosen parameters in the known families. That
is what the score measures, cell by cell, against the smallest certified
known curve of each cell.

ICARM ranks *known* curves by rank against size per torsion group ("can we
find small curves of high rank?"); its curve pages list a regulator, a
verification date and a free-text commentary ("BSD + GRH certified to rank
31" for #302) without describing the verification procedure. Epoch AI's
FrontierMath lists rank ≥ 30 as an
[open problem](https://epoch.ai/frontiermath/open-problems/elliptic-curve-rank)
(`tools/export_epoch.py` writes their JSON). This challenge plays ICARM's
game with certification added: **find your own curve**, prove its points, and
make it smaller than the record of its torsion group and rank. The records
themselves are certified too: a published curve defines a record only if this
verifier certifies its published points.

This is a Yukon optimization benchmark (Eigen / Layr-Labs), a sibling of the
[Heesch challenge](https://github.com/Layr-Labs/heesch-challenge) and built on
the same contract: plain-text submission, everything recomputed from scratch,
fail-closed gates, and a score that can only move with a proof.

---

## The benchmark, precisely

You are given a Python harness (stdlib-only) driving a pinned **PARI/GP
2.17.4** that:

1. **Parses** every `submission/board/*.ec`: one curve per file, five integer
   a-invariants and a list of rational points (format below). At most 512
   files; non-`.ec` files in the directory are ignored.
2. **Verifies** each file — trusting nothing you claim:
   - the discriminant is nonzero and the model is minimal at every prime
     `p ≤ 10⁶` dividing it (`CURVE_SINGULAR`, `CURVE_NON_MINIMAL`);
   - the curve is **new**: its `|a_p|` fingerprint (168 primes) matches none
     of the 5586 curves frozen from ICARM, Dujella's tables and the LMFDB on
     2026-10-15 — which also catches isogenous curves, quadratic twists and
     re-models (`KNOWN_CURVE`);
   - every point satisfies the equation **exactly**, no point is repeated or
     listed together with its negative (`POINT_OFF_CURVE`, `POINT_DUPLICATE`);
   - points of finite order are dropped (exactly, via `ellorder`) — they
     never count;
   - the remaining points are **certified independent**: the canonical height
     pairing matrix is computed at two precisions (57 and 115 digits), the
     entries must agree to `10⁻⁴⁰`, and Weyl's inequality then proves it
     positive definite from its smallest eigenvalue (`POINTS_DEPENDENT`, with
     the exact relation named when one is found). The count of surviving
     points is `certified_rank`.

   Two board files in one class up to isogeny and quadratic twist are
   rejected (`DUPLICATE_CURVE`): by policy the board counts a class once, so
   it holds one curve per isogeny-and-twist class, and every quadratic twist
   of a known curve is ineligible whatever its rank (`RULES.md` §2). Any
   failure rejects the whole submission, with the file name in the detail.
3. **Places** each curve in its **cell**: its torsion group (`elltors`, one
   of Mazur's 15) and its certified rank, keyed like `"2/9"` or `"6x2/4"`.
4. **Scores** the board as

   ```text
   cell          = (torsion group, certified rank)             15 Mazur groups × rank
   contribution  = 100·(baseline − H)/baseline   if H < baseline   (0 otherwise)
                 = 100 + 100/(1 + H/50)          if no known curve of that torsion reaches the rank
   score         = Σ over cells of the contribution of the board's smallest curve there
   ```

   `H = (1/12)·log max(|c4|³, c6²)` is the naive height of the minimal model
   (unchanged by the translations that keep a model minimal; ICARM's
   leaderboard plots `12·H`). The **baseline** of a cell
   is the smallest certified known curve of that torsion whose certified rank
   is **at least** the cell's (a rank-13 curve can set the rank-12 baseline),
   frozen in `ec_verify/records.json` and never changed; cells with no certified known
   curve have none and are breakthroughs, worth more than 100 and more for a
   smaller curve. A certified rank of 0 proves nothing: such a curve is
   accepted and listed but contributes 0 and never holds a cell. The
   submission's own improvement, `metrics.improvement`, is its score minus the
   previous board's total. The result is written to `score.json`.

One linear leaderboard: any new record anywhere promotes
(`minScoreImprovementBips: 0`), because the total rises by exactly the new
record's contribution; a submission that verifies but sets no record is
accepted and not promoted. The board accumulates across submissions, so a
submission that regresses a held cell — removes or enlarges the curve that
holds it — is rejected (`BOARD_REGRESSION`).

### What "valid" means

A submission scores only if all of the following hold for every board file;
otherwise it is rejected with a stable code, exit 1, and no score file:

- **A curve.** Integral a-invariants (≤ 400 digits each), nonzero
  discriminant, minimal at every prime `≤ 10⁶` dividing the discriminant.
- **A new curve.** Not isomorphic, isogenous or twist-equivalent to any curve
  in the frozen known set (`ec_verify/known_curves.json`), nor to another file
  on the board; good reduction at 70 or more of the 168 primes below 1000, so
  that the fingerprint can decide (`CURVE_TOO_MANY_BAD_PRIMES`).
- **Points on it, exactly.** Rational arithmetic, no floats. Duplicates and
  negatives are rejected rather than silently merged.
- **Certified independence.** Torsion is stripped; the height pairing must be
  positive definite by the certificate of `docs/verifier-architecture.md` §3
  (two-precision agreement, Weyl margin). A near-zero regulator is a
  rejection, not a small number.
- **Within caps.** ≤ 64 points, ≤ 2000 digits per coordinate, file ≤ 1 MiB,
  regular file only, ASCII; at most 512 files, named `[A-Za-z0-9._-]+.ec` (no leading dot),
  verified within a 20-minute budget for the whole board.
- **No regression.** Every cell the previous board holds is held by this
  board with a curve at least as small.

An empty board is valid and scores 0. `#CLAIM rank=k` is recorded and
compared, never trusted: claiming 31 and certifying 30 places the curve in
the rank-30 cell and notes the discrepancy.

### Reference numbers

| | rank | curve |
|---|---|---|
| Largest rank known exactly (unconditionally) | 22 | Samuel Li, September 2026 (21: Li; 20: Elkies–Klagsbrun 2020) |
| Class record (lower bound) | ≥ 31 | ICARM #302, Alpöge–Howell with Claude, 2026-08-23 |
| **This benchmark's baseline** | — | the known records (board total 0); the first new record in any cell promotes |
| Where this challenge is played | 1–31 | the 156 cells with a certified known record (every torsion group, ranks 1 to 31 for trivial torsion, 1 to 3 for `Z/8 × Z/2`), plus every cell above them |

### The ladder — what is actually feasible

Measured during the design pass (PARI/GP 2.17.4, one laptop core):

| rung | how | cost |
|---|---|---|
| **5** | 20 000 random small curves sieved by `S(10³) = Σ_{p≤10³} a_p·log p / p`, 2-descent on the best 40 | 10 s |
| **6** | 300 000 curves, two-stage sieve `S(10³)` then `S(10⁵)`, 2-descent on the best 150 | 3 min |
| **8–12** | sieve wider, use `ellrank`: it settles the rank quickly and finds *most* of the generators, not all (on the repo's own rank-14/15/17 witnesses, `tests/corpus/rkeq*.ec`, its points span an index-31 sublattice or miss one generator); saturate with `ellsaturation` and hunt the rest yourself | hours |
| **18** | specialise a published rank-18 family over ℚ(t) (Elkies 2006; Dujella's [families table](https://web.math.pmf.unizg.hr/~duje/tors/generic.html)) and certify its sections | a day of transcription |
| **19–24** | sieve specialisations by `S(X)` and hunt extra points — the road Fermigier (22, 1997) and Martin–McMillen (24, 2000) took | weeks of compute |
| **25–31** | new families (Elkies 2006 → 28; 2024–26 → 29, 30, 31) | research |
| **32** | — | a research result, outside this benchmark's target |

The sums `S(X)` are a search heuristic for your own sieve; the verifier does
not compute them and they do not score. The torsion-group families
(Dujella's tables list one per Mazur group) work the same way with a
prescribed torsion point.

Progress is a smaller curve in any cell. Every cell has its own record, so
the same search effort is better spent where records are thin: a rank and
torsion group with few known curves, or a family parameter that gives a
smaller model at a rank you can certify. A curve certified at a rank no known
curve of its torsion reaches is a breakthrough worth more than 100. Known curves score
nothing, so nobody can paste their way up, and listing fewer points simply
places a curve in a lower cell.

Measured costs (Apple M-series, PARI 2.17.4): certifying the 31-point record
witness takes ~0.8 s, and every gate takes well under a second per curve on
the known record curves, so a board of a few hundred files verifies in minutes; the harness budget for the whole board
is 20 minutes and is the operative limit (the file cap is 512). Per-submission cost is cents.

---

## How to play

1. Install:

   ```bash
   ./setup.sh                # venv + pinned PARI/GP (uses a matching `gp` on PATH, else builds it)
   pip install -e ".[test]"  # optional: run the tests
   ```

2. Write each curve and its points into its own file on the board:

   ```bash
   $EDITOR submission/board/my-curve.ec
   ```

   `submission/board/README.md` restates the format rules; it is not a `.ec`
   file and is ignored. The curve must be your own find: see
   [`RULES.md`](RULES.md) before you add anything you did not search for
   yourself.

3. Verify locally:

   ```bash
   python -m ec_verify submission/board/my-curve.ec   # one curve: its cell, the record, the percentage
   python -m ec_verify submission/board               # the board as the harness scores it
   python -m ec_verify submission/board --parent /path/to/previous/board   # with the regression gate
   ```

   Or import it: `from ec_verify import verify_text` returns a `CurveResult`
   (its `cell` and `naive_height`) or raises `VerifyError` with a stable
   `.code`. In a search loop, `--no-novelty` skips the known-curve gate;
   never submit without it.

4. Submit via Yukon. The benchmark job re-runs everything in a sandbox,
   checks the board against the previous frontier, and writes

   ```json
   { "score": 0.0, "metrics": { "board_total": 0.0, "parent_total": 0.0, "improvement": 0.0, "scoring_version": 6, "records_version": 5, "parent_gate": "checked", "curves": [], "cells": { "1/1": { "baseline": 0.939132, "board_best": null, "contribution": 0.0, ... }, ... }, "improvements": [], ... } }
   ```

   That is the empty board. Each accepted file appears under `curves` with
   its `cell`, `certified_rank`, `naive_height`, `contribution` and
   `holds_cell`; each cell the board gained appears under `improvements` with
   its `delta` in points.

The participant walkthrough (format, every error code, strategy) is
[`docs/submitting.md`](docs/submitting.md); the working contract for coding
agents is [`AGENTS.md`](AGENTS.md).

### File format

```
E <a1> <a2> <a3> <a4> <a6>       # integers: y² + a1·xy + a3·y = x³ + a2·x² + a4·x + a6
P <x> <y>                        # one rational point per line, num/den form (den omitted ⇒ 1)
P 58446054919170749975942104376446/9 -289145377197241504032247540119122580900747897469/27
#CLAIM rank=30                   # optional; recorded, never trusted
# comments and blank lines are ignored
```

One curve per file, directly under `submission/board/`, named
`[A-Za-z0-9._-]+.ec`. Integer tokens are ASCII `-?[0-9]+` only.

### What you can edit

You may modify **anything inside `submission/`**; the scored files are
`submission/board/*.ec`, and anything else under `submission/` (your search
tooling, notes) is ignored by the verifier. You may **not** touch the harness
— `ec_verify/`, `harness/`, `benchmark.sh`, `benchmark.json`, `setup.sh`,
the workflows — which re-verifies everything from scratch; `editablePaths` in
`benchmark.json` is the source of truth and the backend re-enforces it
server-side. Everything under `submission/` becomes public if promoted, and
the promoted board is the next submitter's starting point: keep every file
that holds a cell. Record progress often with `yukon notes add` (public —
strip secrets).

### Rules of the competition

The verifier enforces what it can; the terms in [`RULES.md`](RULES.md) cover
the rest. In short: every curve you submit must be one you found. Copying a
curve from the ICARM leaderboard (a similar competition, where curves belong
to the participants who submitted them), from Dujella's tables, the LMFDB,
the literature, or from another participant here, and submitting it to this
challenge is prohibited — including as an isogeny, twist or re-model, and
regardless of whether the frozen `KNOWN_CURVE` gate catches it. Specialising
a published *family* at your own parameter is allowed. After the competition
closes the maintainers audit every promoted curve against those sources and
against submission histories; a copied curve is struck from the board and any
reward distributed for it is revoked.

---

## Layout

| Path | What it is |
|---|---|
| `submission/board/*.ec` | **The scored files** (anything under `submission/` is editable) |
| `ec_verify/` | The verifier: parser, exact point checks, novelty fingerprint, torsion, certified regulator |
| `ec_verify/cells.py` | Torsion and cell keys, the contribution formula, the board caps |
| `ec_verify/board.py` | The board: per-file verification, duplicates, cells, the regression gate, improvements, the total |
| `ec_verify/records.py`, `ec_verify/records.json` | The frozen cell table: record and baseline of every (torsion, rank) cell |
| `ec_verify/known_curves.json` | The frozen known-curve ledger: fingerprints, torsion, heights, certification status |
| `harness/verify.py` | The Yukon evaluator: grades `submission/board/`, writes `score.json` |
| `benchmark.sh`, `setup.sh`, `benchmark.json` | Sandbox and parent-board export, environment (pinned PARI/GP), manifest |
| `tools/fetch_known.py`, `tools/certify_known.py`, `tools/build_known.py`, `tools/build_records.py` | The record-store pipeline: fetch the sources, certify their points, build the ledger, build the cell table |
| `tools/fetch_records.py`, `tools/export_epoch.py` | Parse Dujella's record pages; export Epoch's JSON |
| `tests/` | Grammar, gates, certificate, cells, board, record store, harness e2e; `tests/corpus/` holds the rank 20–31 record witnesses |
| `RULES.md` | The terms of the competition: originality, what counts as copying, the post-close audit |
| `docs/submitting.md` | Participant guide |
| `docs/verifier-architecture.md` | Normative: pipeline, the certificate and its one assumption, cells and the board total, budgets |
| `docs/records.md` | The record store: files, schemas, the freeze pipeline (for maintainers) |
| `docs/THREAT-MODEL.md` | Security model: assets, attacker, controls, probed attacks, residual risk |

## Credits

Record curves and their points are from Andrej Dujella's
[tables of high-rank elliptic curves](https://web.math.pmf.unizg.hr/~duje/tors/tors.html)
and the [ICARM Elliptic Curve Rank Leaderboard](https://elliptic-rank.icarm.cloud/)
(4235 curves, the whole database, frozen 2026-10-15 as the known set and the
source of most records), every curve on Dujella's rank and per-torsion record
pages (1257), and 94 curves from the [LMFDB](https://www.lmfdb.org/): eight classics
and the smallest curves of every torsion group and rank in Cremona's tables. All arithmetic is PARI/GP (the PARI Group,
Université de Bordeaux), pinned by version and tarball SHA-256.
