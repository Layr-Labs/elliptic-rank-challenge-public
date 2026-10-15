# Submitting: the participant guide

Everything you need to go from a curve to a scored submission. The verifier
is deterministic and fail-closed: the same files always score the same, and
nothing scores without a certified witness. Normative details live in
`verifier-architecture.md`; this page is the practical path.

## 1. The loop

```bash
./setup.sh
python -m ec_verify path/to/curve.ec            # one curve: cell, record, percentage (< 1 s)
python -m ec_verify submission/board            # the board as the harness scores it
python -m ec_verify submission/board --parent /path/to/previous/board   # with the regression gate
python -m pytest tests -q
```

1. Write each curve in its own file under `submission/board/` (§2).
2. Iterate on single files, or in code with `verify_text(text, VerifyConfig(conductor_enabled=False))`
   — same gates, same certificate. The result carries `cell` and `naive_height`; the
   single-file CLI adds the cell's `cell_record` and your `contribution`. In a tight search loop `--no-novelty` (or
   `novelty_enabled=False`) skips the known-curve gate; run the full check before submitting.
3. Run the board once, ideally with `--parent` pointing at a checkout of the current
   frontier's `submission/board`, so the regression gate runs as it will in CI.
4. Submit via Yukon. The benchmark re-verifies every file from scratch in a sandbox, gates
   the board against the previous frontier, and writes `score.json`; any new record in any
   cell raises the total and promotes.

## 2. The file

One curve per file, directly under `submission/board/`, named `[A-Za-z0-9._-]+.ec` with no
leading dot. At most
512 files. Files that do not end in `.ec` (such as the directory's `README.md`) and
subdirectories are ignored; a `.ec` entry that is not a regular file, or has a name outside
that pattern, is a rejection.

```
# comments and blank lines are ignored
E 1 0 0 -201769035260418549083594900060734240952308696994802735114305555 1151107939141058565733479426024323225135665982951300586808823640527729578307228357301072889377
P -4761204159891138283979053265906 44764265461782973805868732003346421827415264953
P 58446054919170749975942104376446/9 -289145377197241504032247540119122580900747897469/27
...
#CLAIM rank=30
```

- `E a1 a2 a3 a4 a6`: integers, the model `y² + a1·xy + a3·y = x³ + a2·x² + a4·x + a6`.
  It must be nonsingular and minimal at every prime `≤ 10⁶` dividing the discriminant
  (the verifier checks with Tate's algorithm prime by prime; it never factors Δ).
- `P x y`: one rational point per line, `num/den` (denominator omitted ⇒ 1), in lowest
  terms or not. Order is irrelevant. Torsion points are dropped, not rejected.
- `#CLAIM rank=k`: optional; recorded as `claimed_rank`, compared, never trusted.
- Grammar is ASCII; fields are separated by spaces or tabs; integer tokens are `-?[0-9]+`.
  Caps: ≤ 400 digits per a-invariant, ≤ 2000 digits per coordinate component, ≤ 64 points,
  ≤ 1 MiB, regular file; everything under `submission/` together ≤ 16 MiB (Yukon's submission
  cap covers the whole editable tree, notes and logs included, and is enforced here too so the
  laptop sees the same rejection). `#CLAIM rank=` takes 1–6 digits.

## 3. The score

```text
cell          = (torsion group, certified rank)             15 Mazur groups × rank
contribution  = 100·(baseline − H)/baseline   if H < baseline   (0 otherwise)
              = 100 + 100/(1 + H/50)          if the cell has no baseline (breakthrough)
score         = Σ over cells of the contribution of the board's smallest curve there
H             = (1/12)·log max(|c4|³, c6²)
```

**`certified_rank`** is the number of your points that survive: on the curve exactly, of
infinite order, and jointly independent by a certified positive-definite height pairing.
It is a lower bound on the Mordell–Weil rank with no GRH, BSD or descent assumption; the
certificate and the one numerical input it relies on (PARI's canonical height, cross-checked
at two precisions) are stated in `verifier-architecture.md` §3. With the torsion group from `elltors`
(one of Mazur's 15, keyed `1`, `n` or `nxm`) it places the curve in exactly one cell, keyed
like `2/9` or `6x2/4`. A board curve counts in exactly one cell. Listing fewer points places
it in a lower cell; that is allowed and it counts only there, so a rank-14 curve that beats
the rank-13 record but not the rank-14 one is best submitted with 13 points. The same curve
twice, at two ranks, is `DUPLICATE_CURVE`. A certified rank of 0 proves nothing: the curve is
accepted and listed but contributes 0, never holds a cell, and is ignored by the
regression gate and the improvements report.

**`H`** is the naive height of your minimal model, rounded to 6 decimals (ICARM's
leaderboard plots `12·H`). Smaller is better.

**Baseline and record.** `ec_verify/records.json` gives every cell that has a certified
known curve its `record`: the smallest certified known curve of that torsion that *reaches*
the rank at the freeze of 2026-10-15, meaning certified at that rank or higher (a rank-14
curve listed with 13 of its points certifies rank 13). Its height is the cell's baseline, and
it never changes; curves published later reach the board only as ordinary submissions. A known curve sits in the
cell of the rank its published points certify, not the rank it claims; curves whose points
do not certify define no record. A cell with no certified known curve has no baseline: a
curve there is a breakthrough worth `100 + 100/(1 + H/50)`, above 100 and strictly more for
a smaller curve.

**Improvement.** `metrics.improvement` is your score minus the previous board's total under
the same table, the number Yukon's promotion depends on; `metrics.improvements` lists each
cell where your contribution exceeds the previous board's, with the `delta`. Without
`--parent` both are `null` and every scoring cell is listed. Because the total is a sum of
per-cell contributions against frozen baselines, it rises by exactly the contribution of
each new record, and any record in any cell promotes.

**Holding a cell.** A board curve holds its cell when its contribution is positive: `H`
below the baseline, or any curve in a breakthrough cell (rank 1 and above). The promoted
board is the next submitter's starting point, and the harness rejects a submission in
which any cell held by the previous board is missing or has a larger best curve
(`BOARD_REGRESSION`). Removing a curve that holds nothing is allowed, and so is certifying
more points on a holder: the curve moves to the higher-rank cell, and its old cell may empty
as long as the board total does not fall (checked once for the whole board, so one new cell's
gain cannot pay for two vacated cells).

**Versions.** `score.json` reports `scoring_version` (6: this formula) and `records_version`
(the frozen table). Neither the table nor the known set changes after the freeze, so a board
always sums to the same number.

## 4. Known curves are rejected — the novelty gate

5586 curves were frozen on 2026-10-15 from the ICARM leaderboard (4235, the whole database),
Dujella's rank and per-torsion record pages (1257) and the LMFDB (94: the classic anchors and the
smallest curves of every torsion group and rank in Cremona's tables) (`ec_verify/known_curves.json`). A
submission whose `|a_p|` fingerprint (168 primes `p ≤ 1000`, bad primes skipped) matches any
of them is `KNOWN_CURVE`: the two agree at every prime good for both and at least 30 primes
were compared. Primes dividing either discriminant are skipped, so a submission must itself
have good reduction at 70 or more of the 168 primes below 1000 (`CURVE_TOO_MANY_BAD_PRIMES`
otherwise; the sharpest known case, icarm-51 twisted by its 83 smallest good primes, still
has 74 good primes, and the thinnest known curve has 139). That floor is what makes the gate
sound: a twist by a smooth `d` is bad at every prime dividing `d`, but it is still compared
at ≥ 70 − 29 = 41 primes against any frozen curve (each has ≥ 139 good primes).
Because `a_p` is an isogeny invariant and `|a_p|` a twist
invariant, this also rejects isogenous curves, quadratic twists and any re-model. The same
fingerprint test runs between your own board files (`DUPLICATE_CURVE`). Specialising a
published *family* over ℚ(t) at a new parameter is allowed and is the intended on-ramp.
Every known curve feeds the gate, whether or not its points certify. The set is frozen; the tools that built it are described in `docs/records.md`.

**The gate is a floor, not the rule.** Because the set is frozen, a curve published after
2026-10-15 — on the [ICARM leaderboard](https://elliptic-rank.icarm.cloud/), in a paper, in the
LMFDB — or taken from another participant of this challenge verifies as new. Submitting it is
still prohibited: every curve on your board must be one you found yourself (`RULES.md` §1).
ICARM is a similar competition and its curves belong to the participants who submitted them
there. After the competition closes the maintainers audit promoted curves against those
sources and against submission histories; a copied curve is struck from the board and any
reward already distributed for it is revoked (`RULES.md` §4).

## 5. The ladder — how the competition actually progresses

| rung | how | measured cost |
|---|---|---|
| **5** | 20 000 random short-Weierstrass curves sieved by `S(10³)`, `ellrank` on the top 40 | 10 s |
| **6** | 300 000 curves, two-stage sieve `S(10³)` then `S(10⁵)`, `ellrank` on the top 150 | 3 min |
| **8–12** | wider sieves and families; `ellrank` settles the rank quickly and finds *most* generators, not all (on `tests/corpus/rkeq14.ec`, `rkeq15.ec` and `rkeq17.ec` its points span an index-31 sublattice or miss one generator), so saturate (`ellsaturation`) and hunt the missing points | hours |
| **18** | specialise a published rank-18 family over ℚ(t) (Elkies 2006; see Dujella's families table) and certify its sections | a day |
| **19–24** | sieve specialisations by `S(X)`, then hunt extra points (Fermigier 1997 → 22, Martin–McMillen 2000 → 24 this way) | weeks |
| **25–31** | new families (Elkies 2006 → 28; 2024–2026 → 29, 30, 31) | research |
| **32** | thirty-two certified points; human review before announcement (`verifier-architecture.md` §7) | research; not the target |

`S(X) = Σ_{p≤X} a_p log p / p` is a search heuristic for your own sieve; the verifier does
not compute it. Progress is a smaller curve in any cell. Every (torsion, rank) cell has its
own record, so a better parameter in a known family, at a rank and torsion group you can
certify, is the plausible game; a cell no known curve of that torsion reaches is worth more
than 100, and more for a smaller curve.
A board holds up to 512 files and each verifies in seconds; the 20-minute board budget is
the operative limit.

## 6. Rejection codes

Stable API — parse them in your loop. `REJECTED: <CODE>: <detail>` on stdout, exit 1,
no score file. A failure in any board file rejects the whole submission; the detail starts
with `board/<file>:`.

| code | meaning |
|---|---|
| `SUBMISSION_NOT_REGULAR_FILE` | a board entry that is not a regular file (symlink, FIFO, device), or `submission/board` itself is not a directory |
| `PARSE_SYNTAX`, `PARSE_EMPTY` | malformed file / no `E` line; non-ASCII; a board file name outside `[A-Za-z0-9._-]+.ec` or starting with a dot |
| `CURVE_TOO_LARGE`, `POINT_TOO_LARGE`, `TOO_MANY_POINTS` | digit / count caps |
| `CURVE_SINGULAR` | discriminant zero |
| `CURVE_NON_MINIMAL` | not minimal at a prime `≤ 10⁶`; the message names `p` and the scale factor |
| `KNOWN_CURVE` | matches a frozen known curve (or an isogeny/twist/re-model of one); the message names it |
| `CURVE_TOO_MANY_BAD_PRIMES` | good reduction at fewer than 70 of the 168 primes below 1000; fingerprint matching (the novelty and duplicate gates) cannot decide such a curve, so the floor holds even with `--no-novelty` (§4) |
| `POINT_OFF_CURVE` | a point does not satisfy the equation (exact arithmetic) |
| `POINT_DUPLICATE` | a point is repeated, or listed with its negative |
| `POINTS_DEPENDENT` | the height pairing is not certified positive definite (`verifier-architecture.md` §3); the message names an exact relation when one is found |
| `HEIGHT_PRECISION_DISAGREE` | heights at 57 and 115 digits differ by more than `10⁻⁴⁰`, or `h(2P) ≠ 4h(P)` to that tolerance, or gp ran below the requested precision (would indicate a PARI problem; never seen) |
| `DUPLICATE_CURVE` | two board files lie in one class up to isogeny and quadratic twist, which the board counts as one curve (`RULES.md` §2); the message names both |
| `TOO_MANY_CURVES` | more than 512 `.ec` files on the board |
| `BOARD_REGRESSION` | a cell the previous board holds is missing from this board or held by a larger curve (the message names the cell and the parent file); or a previous-board file no longer verifies at all, which means the environment changed under the frontier and needs a maintainer; or, on the official runner only, the previous board could not be exported at all (the message says so and names `EC_PARENT_REF`): that is a runner fault, not your board, so report it rather than changing files |
| `RESOURCE_EXCEEDED` | file cap (1 MiB) or the 16 MiB cap on everything under `submission/`, a stage timeout (heights 300 s, torsion and minimality 120 s, fingerprint 60 s), or the 20-minute board budget (a stage cut short by it says so) |
| `GP_UNAVAILABLE`, `GP_VERSION_MISMATCH` | environment: the pinned gp is missing or wrong (not your bug; run `./setup.sh`) |
| `GP_ERROR` | gp reported an error while working on the file (the message quotes gp's `***` line, e.g. a PARI stack overflow on a pathological input), or its output could not be read (`internal error: …`); the detail names the file. Reproducible with `python -m ec_verify` on the same file; if it is not, report it |

## 7. What `score.json` tells you

`score` is the board total. Under `metrics`:

- `scoring_version`, `records_version`, `known_frozen_at`, `pari_version`, and
  `parent_gate` (`checked` when the previous board was exported, `skipped` for a local run
  without `--parent`; the official run of a submission is always `checked`, anything else is
  rejected), plus `board_total` (equal to `score`).
- `curves`: one entry per board file, with `file`, `ainvs`, `cell`, `torsion`,
  `torsion_key`, `torsion_order`, `certified_rank`, `naive_height`, `naive_height_source`
  (`12·H`, ICARM's normalisation), `naive_height_log10`, `discriminant_digits`,
  `fingerprint`, `points_submitted`, `points_torsion_dropped`, `regulator`,
  `regulator_lower_bound`, `height_lambda_min`, `height_precision_agreement`,
  `claimed_rank`, `claim_discrepancy`, `verified_claim` (a sentence stating exactly what was
  established), `cell_record` (the cell's known record, or `null`), `contribution` (what
  the curve alone would earn in its cell; 0 at rank 0), `holds_cell`,
  `pari_version`, and `conductor` (informational: a decimal string, attempted only when
  |Δ| has ≤ 40 digits, else `null`).
- `cells`: every cell with a board curve or a known record, with `torsion`, `rank`,
  `known_record`, `baseline`, `board_best` (`file` and `naive_height`, or `null`) and
  `contribution`.
- `improvements`: one entry per cell where the board's contribution exceeds the previous
  board's, with `cell`, `file`, `H`, `contribution`, `parent_contribution` and `delta`; plus
  `parent_total` and `improvement` (`null` without a previous board) (§3).
- `ignored_parent_files`: previous-board entries the loader skipped (bad name, not a regular
  file, oversize), with their codes; the same loader would have rejected that board as a
  submission, so they never held a cell. A previous-board file that loads but no longer
  verifies, for any reason, rejects the run (`BOARD_REGRESSION`).

Everything except `conductor` is byte-for-byte reproducible on any machine with the pinned
gp; the conductor is bounded by a wall-clock cap and never enters the score.

## 8. Troubleshooting

- **`POINTS_DEPENDENT` on points you believe independent:** check the named relation
  with `elladd`/`ellmul`; a hidden dependency at record scale is usually a transcription
  error in one coordinate. Denominators must be exact.
- **`KNOWN_CURVE` on a curve you found yourself:** you rediscovered a published curve (or a
  twist/isogeny of one); the message names it. Move on — the fingerprint is not negotiable.
- **`CURVE_NON_MINIMAL`:** run `ellminimalmodel` (or `elllocalred` at the named prime) and
  transform your points with the returned change of coordinates.
- **`BOARD_REGRESSION`:** you removed or enlarged a curve that holds a cell; restore it.
- **Accepted but not promoted:** no cell's contribution exceeds the previous board's;
  `improvements` is empty. Check `contribution` and `holds_cell` per curve.
- **Same files, different machine, same score?** Yes. If not, report it.
