# The budget pattern (owner, 18 Nov 2026, 20:05 UTC): which curves can be high-rank AND proved

Second pattern found tonight (the first is `helpers/CONICS.md`). It explains the first one, the "prime-rich"
classes, and why records cluster just above the height where a rank first appears; and it says where NOT to look.

## 1. The measurement (work/patterns/tam_records.gp, bsd_sample.gp; frozen lists as calibration only)
For a curve E let **b = Omega * prod(c_p) / #T^2** (real period times Tamagawa product over torsion squared; PARI
`ellbsd(E)`). Birch and Swinnerton-Dyer: L^(r)(1)/r! = b * Reg * #Sha.

log10 b, quantiles 10 / 50 / 90 %:
| family, rank | certified top-rank curves | random coprime pairs of the same height |
|---|---|---|
| Z/8 rank 6 (4 frozen + our 3: -0.23, -1.6, -2.5) | -1.7 / -1.1 / -1.0 | -5.4 / -4.5 / -3.2 |
| Z/8 rank 5 (68) | -2.7 / -1.1 / -0.2 | -5.5 / -3.3 / -2.0 |
| Z/2xZ/6 rank 6 (4) | -1.3 / -1.0 / -0.6 | -5.4 / -4.5 / -3.2 |
| Z/2xZ/6 rank 5 (109) | -2.5 / -1.3 / 0.1 | -5.7 / -3.5 / -2.2 |
| Z/7 rank 6 (1), rank 5 (31) | -1.1; -1.7 / -0.9 / 0.2 | -5.8; -4.8 / -2.8 / -1.5 |
| Z/9 rank 4 (7) | -1.6 / -1.1 / -0.4 | -3.9 / -1.9 / -0.8 |
| Z/10 rank 4 (25); the 2 classes of our complete box 1500: -0.40, +0.47 | -3.4 / -2.4 / -1.0 | -5.0 / -3.0 / -1.8 |
| Z/12 rank 3 (94); Z/2xZ/8 rank 3 (34) | -3.1 / -2.0 / -0.6; -2.4 / -1.6 / -0.2 | -4.9 / -2.6; -3.6 / -2.4 |
In the three K3 families every certified rank-5/6 curve has b between 0.003 and 3, three to four orders of
magnitude above a random curve of its height. b falls with height like 1/h^1.6 (period ~ 13 e^-H, Tamagawa
product grows slowly): the median of log10 b for Z/8 is -3.7, -5.2, -6.8, -7.6 at h = 1e5, 1e6, 1e7, 3e7.

## 2. Why: two effects, both through BSD
1. **Provability.** Reg * Sha = Lstar / b. Calibrated on our own curves where the search found every generator:
   Lstar = 6e4 to 7.5e4 (rank 4, Z/10), 5e4 to 9e4 (rank 5, Z/8), 7e5 (rank 6, Z/8): about 6e4 * 10^(r-4).
   So the regulator, hence the size of the generators, is fixed by b. Our GPU search on 2-coverings reaches
   canonical height 45 / 55 / 64 at B = 1e5 / 1e6 / 1e7, so r generators are findable only if
   b >= b_min(r) = Lstar_r / (0.3 * reach^r):  rank 4: 0.05;  rank 5: 0.011;  rank 6: 2.4e-3;  rank 7: 5.5e-4
   (B = 1e5; a factor 2 to 4 lower at B = 1e6).
2. **Existence.** The chance of extra rank rises with b (the L-value has to vanish, and b is the size of the
   smallest non-zero value it can take): in our Z/8 census below 6e5 the density of rank 6 is 2.6e-8 among pairs
   with b >= 0.08, 5e-10 for b in 0.003..0.08, and under 5e-12 below (0 curves among 2e11).
Both say the same: a new cell must be looked for among b-rich curves, and nowhere else.

## 3. What it explains
- "Prime-rich" classes (Z/8: 7-8 primes in p(p-q)): at a prime of a width-w cusp the reduction is split with
  c_p = w * v (the torsion point generates the component group), so prod(c_p) = 8^(number of such primes).
- Conics: on a split conic the cusp forms factor further into linear forms of the conic's parameter, so conic
  points have more prime factors, hence larger b (first-descent pass rate 2.4 % against 0.05 % in the box).
- Onset clustering: the number of b-rich pairs per octave of height stays roughly constant while the number of
  pairs quadruples, so the records sit at small height (Z/8 rank 6: seven below 6e5, none found from 6e5 to 1e6).

## 4. What it kills (checked tonight)
- Conic sweeps above about 1e8 (b <= 1e-6 for all but smooth parameters): hunt-conic's 7e8 tier is at idle
  priority, the 1e10 plan is dropped. Medians of log10 b of second-descent survivors: -1.9 (1e4-1e5), -3.3, -4.4,
  -5.2 (1e7-1e8, maximum -3.2).
- The two Z/2xZ/6 conic curves with class-wide bound 7: log10 b = -4.4 and -5.9; Magma's deeper descent: rank <= 3.
- The Z/10 candidates with class-wide bound 5 found so far: log10 b from -4.7 to -9.0, i.e. generators of height
  160 to 1100 if the rank were 5. One of them (-94215/77378) has sums like a rank-5 curve; it cannot be proved.
- Z/2xZ/8 rank 4: of 24,219 curves with second-descent bound 4 in the shell 1000-10000 exactly one has
  b >= 0.02, and the 30 best by b have long sums of rank 0-2. Parked for good.
- 12/4 undercut: complete box to 2000, 97 curves below the baseline with class-wide bound 4; best by sums
  (-148/503, H 26.30, one point of height 29.9 found): Magma's deeper descent gives rank <= 2. Closed.
  10/4 undercut: the 53 candidates have sums of rank 2 (max S(3e6) 14.8 against 24-35 for the rank-4 controls). Closed.
- Z/12 rank 5 (onset near H 40): b is below 1e-9 there. Not a target.

## 5. What it enables: budget-first exact sweeps
Enumerate only pairs with b >= b_min(target) and run the exact chain on all of them with NO Mestre-Nagao cut:
they are few (Z/10: about 1e6 per octave of height at log10 b >= -2; Z/8: about 2e9 per half-decade at
log10 b >= -1.9) and they are the only curves that can score. This reaches heights no box sweep can.
- Z/10 rank 5: job E (box 4e4, best score first) finds every provable one below 4e4 tonight; then the same chain
  on b-rich pairs only from 4e4 to 3e5 (hunt-selmer, job F).
- Z/8 and Z/2xZ/6 rank 7: the queued class sweep (n1 >= 9, 4e6-1.4e7) IS the budget-first sweep there
  (b >= 5e-4 needs n1 >= 9-10 at these heights); stays in the GPU queue behind the current shell.
- Late candidates: `tools/triage.gp` (budget, long sums, predicted generator height) before any point search or
  Magma job. Search points on the member of the isogeny class with the largest b (the regulator is smallest there).

## 6. Corrections after the independent review (review/budget/REPORT.md, 21:22 UTC)
- No kill of section 4 is reopened. The correlation of the class-wide bound with b is real on our complete boxes
  (not selected by finding points): Z/2xZ/6 box 2e4, the fraction whose class bound stays at 5 rises from 0.23
  (log10 b < -2.5) to 1.00 (>= -0.5); Z/10 box 1500 from 0.07 to 1.00. The zero density quoted in 2.2 for the bottom
  bin is a provability artefact, not an existence measurement.
- Lstar: use 4e4 (rank 4), 1e5 (rank 5), 1e6 (rank 6), measured on saturated bases; the formula 6e4 * 10^(r-4) was
  ten times too high at ranks 5 and 6. Lattices are not round: the tallest generator is about 1.6 times the
  geometric mean. Thresholds to prune on: rank 5 (Z/10): below 1e-3 (between 1e-3 and 1e-2 search at B = 1e6);
  rank 6 (Z/7): below 1e-3; rank 7: do not prune on b unless log10 b < -5; rank 4: a round lattice is provable
  near 0.01.
- Take b as the maximum over the isogeny class and search that member (Z/10: the fibre is not the best member for
  four curves in five; factors 4, 16, 64). `tools/triage.gp` now prints the class maximum.
- "PARI found no full basis" is not "unprovable": 34 undecided curves of the Z/10 box 1500 have log10 b >= -1.3.
  They go through the all-classes search (proof-covers).
