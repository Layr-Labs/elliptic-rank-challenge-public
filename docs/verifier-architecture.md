# Verifier architecture (normative)

## 1. Acceptance rule

A submission is the set of files `submission/board/*.ec`. It scores iff every file passes
every gate below and the board passes the board gates (§4). Each accepted curve has a
`certified_rank`, established by exact arithmetic plus the certificate of §3, and a torsion
group; together they place it in one cell. The score is the board total: the sum over cells
of the percentage by which the board's smallest curve in the cell undercuts the cell's
frozen baseline (§4). No GRH, BSD, or descent assumption enters any part of it: the torsion
group is exact, the size is computed from the a-invariants, and the rank is a lower bound
whose one numerical input is PARI's canonical height; §3 states what the certificate relies
on and why the margins make that safe.

## 2. Pipeline (cheapest gate first; every failure is a stable code)

| stage | module | computation | trust basis |
|---|---|---|---|
| 0 | `harness/verify.py` | board directory listing, name pattern; per file `O_NOFOLLOW`, `S_ISREG`, bounded read, strict decode | file discipline |
| 1 | `ec_verify/parse.py` | ASCII grammar, caps | our parser |
| 2 | `ec_verify/curve.py` | Δ from a-invariants, Δ ≠ 0; `c4`, `c6` and `H` | exact integers |
| 3 | `ec_verify/points.py` | equation check, duplicate/negative | exact `Fraction` |
| 4 | `ec_verify/curve.py` | trial-divide Δ by p ≤ 10⁶; `elllocalred(E,p)[3][1] == 1` | Tate's algorithm per prime |
| 4b | `ec_verify/fingerprint.py`, `known_curves.py` | `\|a_p\|` for the 168 primes p ≤ 1000 vs the frozen known set; candidate must be good at ≥ 70 of them | isogeny/twist invariance of `\|a_p\|` |
| 5 | `ec_verify/torsion.py` | `elltors` (the torsion group, hence the cell), `ellorder` per point (0 ⇒ infinite) | exact group arithmetic |
| 6 | `ec_verify/regulator.py` | certificate (§3) | numerical: two-precision agreement plus eigenvalue margin (§3) |
| 9 | `ec_verify/curve.py` | `ellglobalred` if \|Δ\| ≤ 40 digits (20 s safety cap); informational, decimal string | never scored |
| 10 | `ec_verify/board.py` | duplicates across files, cells, parent gate, total (§4) | fingerprints; exact comparison of 6-decimal `H` |

Stages 1–9 run once per file (`ec_verify/verify.py`); stage 10 runs once per board. The
numbering keeps the earlier stage numbers; stages 7 and 8 were retired with the old formula
(History).

`ec_verify/gp.py` is the only PARI entry point: scripts are assembled from
validated integers, run as `gp -q -f` with `parisize = 2²⁸`, `parisizemax = 2³¹` and `recover = 0`
passed on the command line (an in-script stack change makes a source-built gp abandon the script),
explicit `realprecision`, a completion sentinel, and per-stage timeouts. With `recover = 0` any
error aborts gp without the sentinel, so the stage fails closed structurally; in addition any
non-warning `***` diagnostic on stderr fails the stage (`GP_ERROR`) naming the error. Warnings
(stack growth, `warning()`), with or without the `*** <function>:` prefix gp adds inside a
function, are not errors. gp's reported version
must equal the pin (`2.17.4`) or the harness rejects with `GP_VERSION_MISMATCH`.

Minimality beyond 10⁶ is not checked: it would require factoring Δ (166 digits at rank
28), and it does not matter — canonical heights are model-invariant, and a model that is
non-minimal at a large prime only has a larger `H`, which can only lower its score. The gate
exists to keep submitted models canonical and small.

## 3. The independence certificate

Let `M = (⟨P_i, P_j⟩)` be the Néron–Tate height pairing matrix of the r points of
infinite order (PARI normalisation). The points are ℤ-independent modulo torsion iff
`M` is positive definite. PARI computes `M` to a requested working precision without a
rigorous error bound; the certificate is:

1. Compute `M₁` at 57 digits and `M₂` at 115 digits. Require
   `max|M₁ − M₂| ≤ ε = 10⁻⁴⁰` (else `HEIGHT_PRECISION_DISAGREE`). Measured agreement on
   every rank 20–30 record curve is `4–6 × 10⁻⁵⁶`.
2. Take `ε` as the per-entry error of `M₂`. Weyl's inequality gives
   `|λ_i(M) − λ_i(M₂)| ≤ ‖M − M₂‖₂ ≤ r·ε`. Require
   `λ_min(M₂) − r·ε > 10⁻³⁰` (`qfjacobi`). Then `M` is positive definite.
3. Report `regulator = det M₂` and `regulator_lower_bound = Π max(λ_i(M₂) − r·ε, 0)`.
4. Two tripwires independent of the matrix: gp must report the requested working precision
   (a lost `default(realprecision)` silently degrades every height), and for every point
   `|h(2P) − 4·h(P)| ≤ ε·max(1, h(P))` — the canonical height is quadratic, so any accuracy
   loss on large points shows up here (measured 10⁻¹¹³ at 115 digits on the record curve).

**The one assumption:** `ellheight` at 115 digits is accurate to far better than
`10⁻⁴⁰`. PARI gives no rigorous error bound for it, so the certificate is not a formal
proof in the sense of interval arithmetic; it is a numerical certificate with a very large
margin. Step 1 is the tripwire: an accuracy failure large enough to matter would have to be
identical at both precisions. The margin is what makes this safe in practice: on the corpus
curves in `tests/corpus/` `λ_min(M₂)` lies between `8·10⁻⁴` (the rank-26 record) and `0.63`
(Dujella's rank-15 curve), and is 0.21 on the rank-31 curve, against the error allowance `r·10⁻⁴⁰` plus the
`10⁻³⁰` acceptance margin of step 2 (which also covers `qfjacobi`'s own rounding), so only a
wrong height, not an imprecise one, could flip the verdict. An exact alternative exists: fix
one prime `ℓ ∤ #E(ℚ)_tors`, choose primes `p` with `ℓ ∥ #E(𝔽_p)`, and identify
`E(𝔽_p)[ℓ] ≅ 𝔽_ℓ` by a chosen generator; the maps `P ↦ (#E(𝔽_p)/ℓ)·P̄` are `𝔽_ℓ`-linear and
kill torsion, so an `r × k` matrix of their values of rank `r` over `𝔽_ℓ` proves independence
with integer arithmetic only. The height certificate is kept because it also yields the
regulator and its margin has been ample on every corpus curve. Dependent inputs behave as expected: the study's own
dependent triple on 5077a gives `λ_min ≈ −10⁻¹¹⁷`, and the verifier names the exact
relation `P₁ + P₂ + 2P₃ = O` after checking it with `elladd`/`ellmul`.

When the certificate fails, the near-null eigenvector is scaled to small integers and
the relation is verified exactly before being reported; an unverifiable hint is omitted.

## 3b. The novelty gate

The known set (`ec_verify/known_curves.json`, 5586 curves frozen 2026-10-15 from the whole ICARM
database (`/database.json`, 4235 curves), every curve on Dujella's `rk*`/`rkeq*` and per-torsion
`z*` pages (1257) and the LMFDB (94: the classic anchors and the smallest curves per torsion and rank
from Cremona's tables))
stores each curve's fingerprint: `|a_p|` for the 168 primes `p ≤ 1000`, `−1` at primes
dividing its Δ. A submission matches a known curve when the two agree at every prime good for both and at
least 30 primes were compared. Soundness: isomorphic and isogenous curves share every `a_p`;
quadratic twists share every `|a_p|`; two curves that are neither agree on 30 values of
`|a_p|` with probability ~∏ 1/(4√p) < 10⁻³⁰, so a false `KNOWN_CURVE` does not occur in
practice. Completeness: every frozen curve has ≥ 139 good primes below 1000 and matches its
own fingerprint (`tests/test_novelty.py`). A twist by a smooth `d` is bad at every `p | d`,
which removes those primes from the comparison; the verifier therefore requires the candidate
to be good at ≥ 70 of the 168 primes (`CURVE_TOO_MANY_BAD_PRIMES`), so at most 168 − 139 = 29
of them are bad for the known curve and ≥ 41 are compared. The invariant
`known_good ≥ 168 − 70 + 30` is asserted over the whole freeze at test time. The floor is a
backstop, not a bound we expect to meet: the sharpest known case, icarm-51 twisted by its 83
smallest good primes (the regression test), still has 74 good primes; a curve found by search
has nowhere near 98 bad primes below 1000.

History: the 2026-09 freeze used `p ≤ 200` (46 primes) and a fixed 30; high-rank curves often
have 17+ bad primes below 200, so four frozen curves (icarm-503, -319, -45, -374) could not be
matched and were accepted with their published points. The 2026-10-05 patch moved to `p ≤ 1000`
and added "compared ≥ half the known curve's good primes"; being an extra condition on matching,
that widened the twist gap it meant to close (icarm-51 twisted by its 83 smallest good primes,
396-digit `a4`, was accepted) and was replaced by the candidate floor above. False negatives are
exactly the curves absent from the freeze; the set is frozen, and curves published later reach
the board only as submissions (`docs/records.md`). Every ledger entry feeds the gate, whatever
its certification status. The ledger also stores each curve's published `rank_claimed`, which is informational: the cell of a
known curve is set by the rank its published points certify (§4). Twists of higher degree (quartic for j = 1728, sextic for j = 0)
change `|a_p|` and are not caught — correctly so: such a twist is a different curve with its own
rank, and its points must be certified afresh.

## 4. Cells, records and the board total

**Cells.** `elltors(E)[2]` gives the torsion structure in PARI's order (`[]`, `[n]` or
`[n, m]` with `m | n`), keyed `"1"`, `"n"` or `"nxm"`; anything outside Mazur's 15 groups is an
internal error. A curve's cell is `"<torsion>/<certified rank>"` (`ec_verify/cells.py`).
`H = (1/12)·log max(|c4|³, c6²)`, rounded to 6 decimals, is its size.

**Records and baselines.** `ec_verify/records.json` (`docs/records.md`) gives every cell with
a certified known curve a `record`, the smallest certified known curve of that torsion with
certified rank **at least** the cell's rank at the freeze (`tools/build_records.py`; the `3/12`
record is a rank-13 curve); its height is the cell's baseline, the reference for the total, and
it never changes. Cells absent from the table are breakthrough cells. A known curve of
certified rank `r` is a candidate record for every cell of its torsion with rank `≤ r`; its
*claimed* rank plays no part, and curves whose points do not certify define no record.

**The total.** With `b` the baseline of a cell and `H` the smallest board curve in it:

```text
c = 100 · (b − H) / b              if the cell has a baseline and H < b
c = 0                              if the cell has a baseline and H ≥ b
c = 100 + 100 / (1 + H / H_REF)    if the cell has no baseline (breakthrough)
H_REF = 50

board_total = Σ over cells (sorted by key) of c
```

`H_REF = 50` sets the slope of the breakthrough term: it is above 100 for every `H` and
strictly decreasing with no clamp (a curve at `H = 50` earns 150, at `H = 10` about 183), so a
smaller curve in a breakthrough cell always raises the total. Cells with no
board curve contribute 0. A certified rank of 0 proves nothing: rank-0 curves are accepted
and listed, contribute 0, never hold a cell, and are ignored by the regression gate and the
improvements report; rank-0 cells are never built (`MIN_SCORED_RANK = 1`). This is
**scoring version 6** (`SCORING_VERSION` in `ec_verify/score.py`); `score.json` reports it
with the `records.json` version. The total is a pure function of the board files and the
table version, strictly increasing under any record improvement in any cell, and is written
as a JSON float without further rounding.

**Board gates.** At most `MAX_BOARD_FILES = 512` files (`TOO_MANY_CURVES`); no two files with
matching fingerprints (`DUPLICATE_CURVE`, so isogenous curves and twists count once); and the
regression gate. Yukon grades a single-parent commit on the previous frontier (dispatched on a
`submissions/<id>` ref; baselines run on `yukon/baseline/*` and have no frontier) and
`benchmark.yml` checks out with `fetch-depth: 2`; on `submissions/*` refs it sets
`EC_PARENT_REF` to the parent commit together with `EC_REQUIRE_PARENT`, which makes the harness
reject a submission run that reaches it without a parent board instead of scoring it with
`parent_gate: skipped`; `benchmark.sh`, outside the sandbox, exports that
commit's `submission/board` blob by blob (`tools/export_parent_board.sh`, so no participant
`.gitattributes` applies) and passes it in as `EC_PARENT_BOARD_DIR`; a set path that is not a
directory is rejected rather than read as an empty parent.
The parent board is verified first with the candidate's settings: a parent file
byte-identical in content to a candidate file (under any name) is verified once and shared;
any other parent file is verified for the gate only, conductor off. A promoted file verified
when its frontier was accepted and, with a frozen table and ledger and a pinned gp, keeps
verifying; one that fails now, whether a gp failure (`RESOURCE_EXCEEDED`, `GP_*`) or a curve or
point code (including `KNOWN_CURVE` after a ledger rebuild), means the environment changed
under the frontier and rejects the run (`BOARD_REGRESSION`), since ignoring it would silently
unhold the cell. Entries the loader skipped (bad names, non-regular files) are listed in
`ignored_parent_files`; the same loader would have rejected that board as a submission.
For every cell the parent holds (a curve with `H` below the
baseline, or any curve in a breakthrough cell), the candidate's smallest curve there must be
no larger (`BOARD_REGRESSION`), unless the parent's holder (matching fingerprint, so the
curve up to isogeny and twist) sits on the candidate board certified at a higher rank. Every
cell vacated that way is then checked once, globally: the candidate's total must not fall below
the parent's, so one destination cell's gain cannot pay for two vacated cells, and whether the
moved curve or a new smaller one is best in its destination does not matter. Without
`EC_PARENT_REF` the gate is skipped and `metrics.parent_gate` says so. benchmark.yml sets
`EC_PARENT_REF` and `EC_REQUIRE_PARENT` together on every `submissions/*` ref, and with the latter
set the harness rejects a run that reaches it without a parent board (`BOARD_REGRESSION`), so no
submission can score ungated. The gate is a gate, not a score input: laptop and CI agree on the
number.

**Improvement.** `metrics.parent_total` is the parent board's total under the same table and
`metrics.improvement` is `board_total − parent_total`; `metrics.improvements` lists each cell
(rank ≥ 1) where the candidate's contribution exceeds the parent's, with the `delta`. This is
the submission's own result; it is reported, not summed.

## 5. Budgets and caps

| item | value |
|---|---|
| a-invariant digits / coordinate digits / points / file / submission tree | 400 / 2000 / 64 / 1 MiB / 16 MiB (`MAX_BOARD_BYTES`, equal to `benchmark.json` `maxSubmissionBytes`; counted over every regular file under `submission/`, as the platform counts it, by `files.check_submission_tree` from `board.verify_board_dirs`) |
| board files / board budget | 512 / 1200 s wall clock for the whole board (`WALL_BUDGET_S`, `ec_verify/board.py`), a hard limit: the harness passes it as `VerifyConfig.deadline` and `gp.run_gp` clips every stage timeout to the time left, so one slow curve cannot carry the run past it; the budget is the operative limit, the file cap sits above the 156 known cells plus breakthroughs |
| PARI stack | `parisizemax = 2³¹` |
| stage timeouts | localred 120 s, fingerprint 60 s, torsion 120 s, heights 300 s, relation hint 20 s (one attempt, ≤ 16 points of ≤ 300 digits), conductor 20 s |
| job timeout | 90 minutes (`benchmark.yml`), the backstop |
| measured (64 dependent points of 1500 digits) | torsion < 1 s, certificate 7.6 s |
| measured (rank-30 record, M-series) | heights 0.2 s |

## 6. Determinism

Same bytes ⇒ same `score` and the same `metrics` on any host with the pinned gp, except
`metrics.curves[*].conductor`, the only machine-dependent field: exact arithmetic for every
gate, fixed `realprecision`, `H` rounded to 6 decimals, cells summed in sorted key order,
sorted JSON keys, no timings in the output, `python -E -s` (no `PYTHON*` variables, no user
site; the stage runs the checked-out tree, the same code as the CLI). Nothing in the output
depends on hash order, so no hash seed is pinned (`-E` would discard `PYTHONHASHSEED` anyway). The conductor is
gated deterministically by discriminant size but bounded by a wall-clock cap, so it may be
`null` on a slow machine; it never enters the score. The harness writes the score into
scratch and `benchmark.sh` copies it out only on success.

## 7. Record procedure

A submission with `certified_rank ≥ 31` is a research result. Before announcement:
re-run the harness on a second machine, re-check the points independently (e.g. Sage or
Magma heights), and confirm the curve is not isomorphic to a known record curve. The
benchmark score stands regardless; the announcement is the human step.


## History

Scoring versions 1 to 4 (2026-09-09 to 2026-09-28) scored a single file as its certified rank
plus a fractional curve bonus: a Frobenius fraction from the Mestre–Nagao sum
`S(X) = Σ_{p≤X} a_p log p / p` at `X = 10⁸`, placed in rank bands frozen in the calibration
table, plus a size fraction measuring the naive height against the median known curve of the
rank (one point at the size record from v3; capped at ½ where fewer than ten known curves
defined the rank in v4). It was replaced on 2026-10-08: the
leaderboard median is not a typical curve at low rank, the band construction scored the
record inconsistently across ranks, and a heuristic for unproven rank has no role in a
smallness contest and favours large curves. The calibration table, the 80-second sum, the
2-descent extra and the single-file layout were removed; per-torsion records against the
smallest certified known curve replaced them.
