# z3sel: fast rigorous rank upper bound for the Z/3 Tate family (3-isogeny descent)

Family: E: y^2 + a x y + b y = x^3, T = (0,0) of order 3. Sieve pair (P, Q) = (b, a^3).
Bound B3 = s_phi + s_phihat - 1 >= rank E(Q) (phi: E -> E/<T>).

## Files
- `search/selmer/z3sel.gp`: `z3bound(P, Q)` returns `[bound, s_phi, s_phihat]`. `z3boundab(a, b)` and `z3loc(a, b)` give per-prime local data.
- `work/z3sel_filter.sh TARGET IN [OUT]`: keeps the `z3` lines with B3 >= TARGET and appends ` B3=<bound>`.
  - The input is copied once, and each bound is joined back by line number with a P,Q check, so a live file that is rewritten mid-run cannot misattach bounds.
  - Appending ` B3=` moves the region tag from `$NF` to `$(NF-1)`.
  - Lines whose first field is not `z3` pass through with ` B3=na`.
  - A z3 line whose bound fails is kept with ` B3=err` (fail-safe). There are 0 such lines in the current rank files.
  - It runs one gp, pinned to core 19 with a 2 GB cap. Stats go to stderr.
- `work/z3sel/`: validation scripts and logs.
  - `val_z3.gp`: stratified random check against the oracle.
  - `t_ranked.gp`, `t2.gp`: ranked curves.
  - `t_pending.gp`, `t_speed.gp`: pending lines and timing.
  - `study3*.gp`: 3-adic study.
  - Snapshots of the input files used here.

Nothing under `search/lottery/`, `queue.tsv` or rankcuts was touched.

## Method (closed form; no point sampling)
The setup is the same as `oracle3.gp`. The Kummer map of the dual isogeny (kernel mu3) is R -> y(R) in Q^*/Q^*3, since the tangent at T is y = 0, and kappa(T) = 1/b.

The formula works on the reduced model: (a, b) ~ (u a, u^3 b), divided out at every l with l | a and l^3 | b. Set D = a^3 - 27 b and S = {3} + primes(b) + primes(D). Every factor above 2^64 is proven prime with isprime.

Local image W_l in Q_l^*/Q_l^*3:

| prime | reduction | W_l | condition on d | h - dim W - 1 |
|---|---|---|---|---|
| l \| b, l ∤ a | split I_{3k} (Tate, q a cube) | everything | none | -1 |
| l \| D, l ∤ b | I_m, T in identity comp. | trivial | d an l-adic cube | 0 (l=2 mod 3), +1 (l=1 mod 3) |
| l \| a, l \| b, l != 3 | additive IV/IV*; E = E_0 + <T>, E_0 3-divisible | <b> | d in <b> (no condition if l = 2 mod 3) | -1 (l=2), 0 (l=1) |
| 3, 3 ∤ a, 3 ∤ b | good ordinary, flat mu3 | units | v3(d) = 0 mod 3 | 0 |
| 3, 3 ∤ a, 3 \| b | split mult. | everything | none | -1 |
| 3, v3(a)=1, v3(b)=0 | | <b>; but units if s = b/(a/3)^3 = 19 mod 27 | | |
| 3, v3(a)=1, v3(b)=1 | | everything if a/3 = b/3 mod 3, else <b> | | |
| 3, v3(a)=1, v3(b)=2 | | everything | | |
| 3, v3(a)>=2, v3(b)=0 / 1 / 2 | | <b> / everything / <b> | | |

At 3, <b> for a 3-adic unit b means: trivial if b = ±1 mod 9 (b is a cube), otherwise units.

- s_phihat = #S - rank_F3(conditions).
- s_phi = s_phihat + sum over l in S of (h_l - dim W_l - 1), by Greenberg-Wiles as in the oracle.
- bound = 2 s_phihat + corr - 1.

Provenance of the rules:
- The rules for l != 3 come from Schaefer's formula plus kappa(T) = 1/b.
- The 3-adic table was read off the oracle on about 7,000 3-adically biased pairs. It depends only on t = b/a^3.
- The good, Tate and flat-mu3 rows agree with theory. The `19 mod 27` exception is empirical, but every check below confirms it.
- Guard: when 3 | a, the rows are used only inside the validated depth, v3(a) <= 25 and v3(a^3 - 27b) <= 40. Beyond that z3bound raises an error, and the filter keeps the line as `B3=err`.
  - Sieve data reaches v3(a) = 9: 3 pending lines.
  - A first guard at v3(a) <= 7 caught exactly those 3 lines. The deep run below then extended the validated range.

## Validation (0 violations, 0 mismatches)
**1. Oracle on ranked curves.** I ran `oracle3.gp` on all 16,098 ranked z3 curves in `logs/lottery/ranked-all.txt` (snapshot 07:00Z; ranks 6-11, r_lo = r_up throughout).
- B3 >= rank on every curve.
- B3 = rank on 15,127 curves (94.0%). The slack B3 - rank is always even:

| slack | curves |
|---|---|
| 0 | 15,127 |
| 2 | 857 |
| 4 | 104 |
| 6 | 10 |

- Tight cases by rank:

| rank | tight / total |
|---|---|
| 6 | 23/23 |
| 7 | 2120/2170 |
| 8 | 10970/11216 |
| 9 | 1766/2145 |
| 10 | 232/504 |
| 11 | 16/40 |

- Named curves (P, Q as in the brief):

| curve | B3 | B3 = rank? |
|---|---|---|
| our 3/10 | 10 | tight |
| our 3/11 | 11 | tight |
| rival 3/10 | 10 | tight |
| rival 3/11 | 11 | tight |

**2. Closed form against the oracle on real data.**
- The closed form equals the oracle on all three numbers for all 16,098 ranked curves.
- It also equals the oracle for all 8,338 pending z3 lines.

**3. Stratified random members: 7,200 in total, plus a deep run of 1,200.** The main runs used two seeds, 9 strata of 300-500 each:
- generic
- smooth b with 4-15 small prime-power factors (repeats allowed)
- 3-adic grid with v3(a) up to 7 and v3(b) up to 20
- near a^3 = 27b with v3(D) up to 19
- l | a and l | b for l in {2,3,5,...,43}
- non-reduced pairs
- tiny D
- 3 | a and 3 | b
- sieve-like

Both signs of a and b were used. Checks on each curve:
- The triple [bound, s_phi, s_phihat] equals the oracle on every curve.
- For every prime of S, the closed-form W_l equals the oracle's Schaefer target plus its sampled image, in the oracle's own c3cls coordinates. There were 0 local mismatches.
- 3.6M extra Q_l-points were sampled at every prime where W_l is not everything. Every Kummer class fell inside the closed-form W_l. This independently tests the cases where the oracle's target is 0 and it samples nothing.
- (-1)^B3 = root number on all 7,200.
- ellrank was run on 893 small curves: B3 >= r_lo on all.

Coverage of the 3-adic types (both runs together):

| 3-adic type | curves |
|---|---|
| good at 3 | 1838 |
| 3 \| b, 3 ∤ a | 1756 |
| v3(a)=1, v3(b)=0: trivial / units | 896 / 656 |
| ... of which s = 19 mod 27 | 201 |
| ... of which v3(D) >= 9 | 134 |
| v3(a)=1, v3(b)=1: everything / <b> | 178 / 180 |
| v3(a)=1, v3(b)=2 | 329 |
| v3(a)>=2, v3(b)=0: trivial / units | 134 / 264 |
| v3(a)>=2, v3(b)=1 | 479 |
| v3(a)>=2, v3(b)=2 | 490 |

Other coverage:
- 2,626 curves have at least one additive prime l != 3: 1,484 instances with l = 1 mod 3 and 2,601 with l = 2 mod 3, including 1,325 with l = 2.
- 2,750 input pairs were non-reduced.

The oracle raised 0 errors ("stuck" or "target too small").

Deep run (seed 31337, 1,200 curves): v3(a) = 8..25 with v3(b) = 0, 1, 2, plus near a^3 = 27b with v3(D) up to 40, partly non-reduced. Results:
- 0 mismatches of the triple and of the per-prime images.
- 0 extra-point failures over 525k points.
- 0 parity failures.

## Speed
| method | time per curve |
|---|---|
| closed form, sieve data (pending / ranked) | 0.07 / 0.08 ms mean; max 1 / 3 ms |
| ... of which factoring b and D | 0.045 ms |
| closed form, random generic pairs (20-digit non-smooth b) | ~1 ms mean |
| `oracle3.gp`, ranked | 5.1 ms mean |
| `oracle3.gp`, pending | 6.6 ms mean, 14 ms max |
| filter script, 8,338 lines (including gp start-up) | about 0.6 s of CPU in total |

The filter's wall time depends on load: 0.65 s with core 19 idle, but 15-40 s when core 19 was contended or the job was parked by the governor at nice 15.

## Pass rates (fraction with B3 >= target; the rest can be skipped without ranking)
| file (snapshot ~07:00-07:20Z) | region | target | lines | B3 >= target | B3 histogram |
|---|---|---|---|---|---|
| rank/pending.txt | z3u10 | 10 | 8338 | 2509 (30.1%) | 7:32 8:5440 9:357 10:2502 11:7 |
| rank/scored.sorted | z3u10 | 10 | 3918 | 528 (13.5%) | 7:26 8:2808 9:556 10:521 11:7 |
| rank/scored.sorted | z3u9 | 9 | 1241 | 40 (3.2%) | 6:3 7:1147 8:51 9:40 |
| rank/scored.sorted | z3u11 | 11 | 561 | 416 (74.2%) | 8:5 9:110 10:30 11:320 12:40 13:56 |
| rank/z3u11-s97e13.deferred.txt | z3u11 | 11 | 1302 | 814 (62.5%) | 8:2 9:439 10:47 11:795 12:16 13:3 |
| (none) | z3u12 / z3w | 12 | 0 unranked lines in any rank file | - | - |

pending.txt held no z3u9, z3u11, z3u12 or z3w lines at snapshot time, so those regions were measured on the other rank files above.

Hindsight on already-ranked curves (ranked-all, by region tag):

| region | ranked | B3 >= target | ranks reached |
|---|---|---|---|
| z3u10 | 1087 | 53 (4.9%) | the only rank-10 curve passes |
| z3u11 | 626 | 331 (52.9%) | none reached 11 |
| z3w | 223 | 196 (87.9%) | B3 >= 12; all are rank <= 11, and the slack is mostly 2-4 (Sha[3] or Selmer excess) |

Since B3 >= rank always holds, no curve of rank >= target is ever dropped.

## Notes for integration
- Usage: `work/z3sel_filter.sh 10 cands.txt kept.txt`. It reads `$2 = P` and `$3 = Q` (Q must be a perfect cube). The format of the rest of the line is free.
- Parity: (-1)^B3 has matched the root number on every curve tested, and every slack on the ranked curves is even. The filter does not use this; it applies only the proven cut B3 >= target.
- The closed form needs only factor(b) and factor(a^3 - 27b). The rows are local symbols, so the same code ports directly to C or the GPU hook (as b3fast did for z2z6) if the sieve ever needs it in-line.
