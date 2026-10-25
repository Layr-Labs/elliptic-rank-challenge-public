# Positive result ready for Yukon submission
Own new curve first-t2x2-p35213-q16986-r13.ec, p35213/q16986,
minimal H25.818611, 13 independently certified points, torsion2x2.
Fresh reproduction regenerates identical witnesses and GP bounds[13,13].
Full board24curves with all23parentfilesbyte-identical:
score111.9791058985735, parent109.84323118123126,
improvement2.135874717342233, parent_gate checked.
Evidence final-board-verified.json, discovery.json, winner-reproduced.json,
winner-search-record.json, SUBMISSION_NOTE.md (12KiB).
Latest frontier f1c42abe0153e017898256f2d1ad9dbe7b47cd32 unchanged afterfetch.
Material dependency43af8780-369f-484f-b8df-47147017654b.
User modelGPT-6.1 Sol/harnessCodex. User authorizedpositive submission.
NEXT: preparedependencies, submitandmonitoractualpromotion/leaderboard.
No needsearchfurther; intendedfirstsmallpositivegoalnowhaswinner.

# Latest update: full-2-torsion X400a research active
No own positive improvement found or submitted yet. User explicitly wants a
positive improvement. User model attribution GPT-6.1 Sol, harness Codex.
Current inherited board: 23 curves, score 109.84323118123126, frontier
f1c42abe0153e017898256f2d1ad9dbe7b47cd32. Dependency Yukon submission
43af8780-369f-484f-b8df-47147017654b. Its new curve2x2-12-own.ec is
inherited and must remain. Parent research/frontier-109-board; use this for
final board regression verification, never claim inherited gain as ours.

Corrected full6 N150000 finished 45000150000 actual pairs, 233.472 seconds.
3977 new rank checks plus4 known, no gain. Initial triangular full6 scan
coverage claim was invalid; corrected square scan and signed validations
replaced it. Current full6 outputs archived under /tmp with hashes.

X400a family reused from inherited file HEADER (published family attributed
to wgxli; ICARM curve4252 is NOT a candidate). Fresh own p/q only:
 A=q(121q-99p)(26p-23q)(75p-44q)
 B=p(21q-74p)(25p-36q)(p+77q)
 E=[0,A+B,0,A*B,0].
A-B=(p-23q)(p-4q)(25p-11q)(74p-121q).
Fixed common-factor primes 2,3,5,7,11,17,19,23,41 from resultant.
x400a-conductor-simd.cpp: full signed box N10000 screened200010000
pairs in8.96s, 30 validation samples exact scalar sieve matches, height
lower bound<=actual GP height with gap<=log2. Deep sieve18088, rank8143
plus1 inherited duplicate; no gain, best fresh rank12 at H23.649257,
not smaller than current rank12 holder18.004452.

CURRENT running: /tmp/ec-rank-x400a-simd 60000 10
output /tmp/ec-rank-x400a-wide-shortlist.txt; progress
/tmp/ec-rank-x400a-wide-progress.txt, tool session31920.
Expected7200060000 counter pairs, ~340s. Next run
.venv-bench/bin/python submission/research/refine-x400a-wide.py
then search.py --input sieve-x400a-wide.txt --limit3000.
Refiner supplies metadata x400a for independently derived point seeds.

x400a-all-sections.gp symbolically solves square factorizations of
x=c*product(subset A factors)*product(subset B factors), degree<=4.
Found24 polynomial sections; four independent at own17/13 (indices1,2,3,7)
height determinant179.58, all exactly on curve. Generated x400a-family.gp
homogenizes x byq^4,y byq^6; search.py now supports its seeds.
These four do NOT prove generic rank5; inherited header claims rank5.

Full notes below are historical and may be outdated.

# Latest update: corrected full-6x2 search active

No positive improvement yet. User explicitly wants continued positive search.
Full-4x2 SIMD 600k completed 179999700000 pairs in 755.34 seconds.
36098 unique shortlist entries were deep-refined, 7725 retained; new rank
checks 5546, no improvement. Remaining eligible height-boundary scan p>600k,
p bounded analytically, screened 5119514652 extra pairs in 29.07 seconds;
3723 new rank checks, no improvement. Full parameter HEIGHT geometry was
covered for cell 4x2/7, but only heuristic finalists had rank checked; do not
claim global rank/height optimality or exhaustive point certification.

CURRENT active experiment is full6-conductor-simd.cpp, executable
/tmp/ec-rank-full6-conductor-simd 150000 8.
Output full6-corrected-shortlist.txt; progress full6-corrected-progress.txt.
Tool session 18049 launched it. At last check denominator 60000, 87 seconds;
expect about 220 seconds total. Output only complete after FINISHED_GRID.
Expected actual counter count N*(2*N+1)=45000150000.

Family from Dujella-Peral arXiv:1306.0027:
 A=37*q^4-84*p*q^3+102*p^2*q^2-36*p^3*q-3*p^4
 B=32*q*(p-q)^3*(p+q)^3*(3*p-5*q)
 E=[0,A,0,B,0], torsion 6x2.
Six singular linear factors q,p-q,p+q,3p-5q,p-3q,3p-q.
Native range now correctly includes p=-N..N and q=1..N; gcd checked.
Do not confuse with full4 triangular domain p>q.

IMPORTANT debugging history: initial full6 port accidentally kept the previous
family's positive triangular p>q loop, though its printed coverage formula
claimed a larger square. That coverage claim is invalid and was discarded.
Signed-modulo indexing was also hardened before activating negative p.
Individual GP rank checks from those old positive candidates remain valid,
so caching their rigorous ineligibility is safe. User was informed candidly.
Corrected version counts counterPairs directly rather than printing a formula.

Validation: validate-full6-simd.py on small N=500 matches scalar table sums
for 50 candidates (25 negative), exactly zero integer counter error.
15 GP trace comparisons agree within 8.9e-15. full6-simd-validation.json.
validate-full6.py separately checks minimal heights and conductor on 100
random primitive p,q, max errors 3.6e-15 and 1.5e-14.
C++ uses /tmp/ec-rank-full6-prime-tables.txt (AP <=2000),
and /tmp/ec-rank-full6-local-tables.txt (minimal scaling and f2 at 2,
indexed by p,q modulo 128). Generated by corresponding .gp files.

Known rank-6 baseline is only a CONTROL, never a candidate. Solving its
j-fiber gave six parameter representations, smallest parameter height 92641.
Recorded in full6-baseline-control.json. All known controls are excluded by
full novelty gate and must never be submitted. Baseline 6x2/6 H22.598355,
8 ledger entries /4 classes; 6x2/5 H14.780907, 113 classes.

Exact Q-isomorphic parameter symmetries were proved algebraically in
validate-full6-symmetries.gp:
 R(v)=2-v; M(v)=(7-v)/(5-3v).
C4,C6 transform by rational scale 16/(5-3v)^2 for M, unchanged for R.
The six normalized parameter pairs are compared lexicographically by
(max(abs(p),q),q,p); only the minimum is retained. This improves shortlist
diversity and loses no curve whose any representative lies in the box.
Proof GP printed FULL6_SYMMETRIES_EXACT.

NEXT: after corrected SIMD completes, run:
 .venv-bench/bin/python -u submission/research/refine-full6-corrected.py
   > submission/research/full6-corrected-refinement-progress.txt 2>&1
 then:
 .venv-bench/bin/python -u submission/research/search.py
   --input sieve-full6-corrected.txt --limit 10000
   > submission/research/search-full6-corrected-progress.txt 2>&1
Refiner handles seven fields (62 p q H S2000 logN mode), dedups models,
GP S10000, both ranking modes, standard row torsion key 6x2.
GP refiner batches are cached under /tmp/ec-rank-conductor-refine-cache.
Rank checks (known / board-duplicate / rigorous upper bound below target)
are cached in /tmp/ec-rank-completed-curve-checks.json, parent SHA-bound.
Search.py uses effort 0 then effort 4 only for viable upper bounds.

Archive housekeeping: bulky older generated files checksum-moved to
/tmp/ec-rank-search-archive-20261025; 24982 completed checks persisted before
archiving old log. Manifest archived-intermediates.json. Public sources remain.
Submission tree was reduced to 3.1MiB before later new data (may need cleanup again).
Do not remove board holders. Research .ec proof witnesses are not board entries.

Rest of older state follows.

# First positive Yukon submission: working state

User explicitly wants a positive improvement and a promoted leaderboard entry.
Do NOT submit a zero-gain test. Exact attribution supplied by user: GPT-6.1 Sol;
harness Codex. No own scoring curve has yet been found or uploaded.

## Environment / execution

Repo /home/linux/ec-rank in Ubuntu WSL. Windows UNC sandbox setup fails.
Working tools: escalated PowerShell process, cwd C:/Windows, invoking
wsl.exe -d Ubuntu --cd /home/linux/ec-rank --exec /bin/bash.
functions store wslExec contains reusable JS wrapper: eval(load("wslExec")),
then await run(bashScript,maxOutputTokens,initialWaitMs). It base64-encodes the
script so PowerShell / WSL quoting is safe. Normal exec and Node REPL cannot start.

Use .venv-bench/bin/python and tools/pari/bin/gp, pinned PARI 2.17.4.
Yukon /home/linux/.local/bin/yukon. Installed g++ supports OpenMP/native SIMD.
Only submission/ edited; verifier, harness, records and setup unchanged.
Scratch binaries and archived generated data are in /tmp. Git fetch only updated refs.

## Frontier

Original board 21 files, verified score 94.11714917597493.
Current promoted frontier 07e2abcc754fd889214af458aeaeb016a2004705,
submission 5efdc131-d52e-4b2c-9962-c0b912aa9de8,
verified score 104.61175471480624, 22 files.
It adds only 4-9.ec, added to local board by exclusive-create after verifying a
separate frontier copy. No local existing file overwritten.
All current board files are inherited. Preserve them; no own board addition yet.

research/frontier-board is the verified 22-file parent for regression checks.
frontier-verified.json and initial-board.json contain that score and metrics.
initial-board-94.json and parent-board preserve the old 21-file snapshot.
Search.py reads initial-board.json to compare against current contributions.

Benchmark ID 0ca81561-012c-48e5-b325-530be7feaeac, name layr-labs/ec-rank.
https://www.yukon.org/ecrank lists promoted submissions only.
Discussions disabled. Installed CLI has no legacy notes add.
Recheck frontier before actual upload: other solvers can advance it.

## Research completed

Independent universal Tate, cached residue-class prime sieves, and K3 searches
have not yet improved a cell. Candidates often have insufficient rank for size.
Published FAMILY formulas only used, never published individual curves.
Sources:
https://arxiv.org/html/2003.00077 equations (3),(4),(7): torsion-2 K3 family.
Published Elkies-Klagsbrun section 11 equations (5),(6): torsion-4 K3 families,
https://msp.org/obs/2020/4-1/obs-v4-n1-p15-s.pdf .
Standard Kubert families, including fresh rational parameters for torsion 5--12.

Own torsion-3 construction: if a1=v-u-w and a3=u*v*w, then (u*v,u^2*v)
lies on y^2+a1*x*y+a3*y=x^3. Pair six linear factors into quadratics and
match products to obtain sections. Inspired by the shape in inherited
board/3-r12-z3k3.ec, so credit that material predecessor in provenance.
Recovered roots (1,5,11), lambda=309/7, K=-13/49, and generated roots
(5,11,19), lambda=3711/13, K=-1189/2197. Wider list in z3-wider-families.txt.
Model: a1=t^2+lambda, a3=K*product(t^2-root_i^2).
Two test fibers of rank 5 are NOT enough to infer generic rank; checking heights
of actual polynomial sections showed most such families only supplied rank 3.
The two listed above each supply five independent sections.

family-section-certificates.json uses FULL official novelty + independence to
certify five points on each own fiber t=101/103. Research witnesses
family-1-sections.ec and family-2-sections.ec have rank 5, heights 13.636121
and 15.892096, score zero. They are under research/, NOT board. Never present
them as improvements. They establish section independence for the family search.

Slow older torsion-3 search processes 89527 and 92286 were stopped, along with
only their verified GP children, after many 45s timeouts. Their source/logs remain.
No current useful search depends on them.
Torsion-4 K3 run search-z4-progress.txt completed: 46 known, 3615 rank checks,
no improvement, 133.4 sec. Isogeny checks of five potential torsion upgrades also
found no improvement.

## ACTIVE: SIMD universal torsion 4x2 search

Family for primitive positive p>q:
y^2=x*(x+(p-q)^2)*(x+(p+q)^2)
=x^3+2(p^2+q^2)x^2+(p^2-q^2)^2*x.
Normally torsion 4x2; official torsion is always recomputed.
These are our chosen specializations, not copied individual curves.

Active source full4-conductor-simd.cpp, executable
/tmp/ec-rank-full4-conductor-simd, launched as:
  /tmp/ec-rank-full4-conductor-simd 600000 12
Output at end: research/full4-conductor-shortlist.txt.
Progress: research/full4-conductor-progress.txt, every denominator 20000.
Tool exec session 16273 launched it. Check if complete before re-running.
12 threads of 16 available logical CPUs. int16 SIMD block counters, SCALE=128.
AP tables primes <=2000 are /tmp/ec-rank-full4-simd-tables.txt,
generated by full4-simd-tables.gp. 303 primes; sentinel 9999 means bad reduction.
Plain-trace counter bounds [-18177,18116], so no overflow.

Two priority modes per target rank:
mode 0 = S(2000) - 0.4*(bad primes <=2000).
mode 1 = S(2000) - 0.5*log(conductor).
These ONLY guide work, not certificates.
Keeps 12000 per mode per cell, with duplicates between modes.
Targets 4x2/6, height <10.740452, and 4x2/7, height <12.948901.
OUTPUT HAS SEVEN FIELDS:
42 p q H plain_S2000 log_conductor mode

Earlier full4-simd.cpp run N=150000, 8 threads screened 11,249,925,000 grid pairs
in 65.79 sec. Its output was GP-refined by prepare-full4-simd.py:
search-full4-simd-progress.txt FINISHED known 85, board duplicate 2,
rank-tested 5966, found False, 87.3 sec. Do not rerun those candidates needlessly.

Fast exact minimal height for primitive p>q:
c4=16*(p^4+14*p^2*q^2+q^4)
c6=64*(p^2+q^2)*(p^4-34*p^2*q^2+q^4)
Hraw=max(log(c4)/4,log(abs(c6))/6).
If p,q both odd, subtract log 2; if p-q OR p+q divisible 8, subtract another log 2.
Otherwise no subtraction. GP still obtains the actual integral minimal model.
30 random checks: max error 1.8e-15, stored in simd-validation.json.

Exact conductor:
log N=log rad_odd(p)+log rad_odd(q)+log rad_odd(p-q)+log rad_odd(p+q)+f2*log 2.
Opposite parity f2=4. Both odd, let d=v2(p-q)+v2(p+q):
f2=3 for d=3, f2=0 for d=4, f2=1 for d>=5.
30 GP global-reduction checks: max log-conductor error 7.1e-15.
SIMD precomputes odd radical logs and small-prime divisor counts through 2*N.

## NEXT STEPS

1. Monitor active 600k SIMD run; meaningful commentary at least every 60 sec.
2. Build GP refinement for its SEVEN-field output. Deduplicate p,q first.
   Obtain integral minimal model with ellminimalmodel. Compute S(5000) or deeper.
   Retain strongest union using both bad-prime bonus and conductor-adjusted score
   per target rank, without accidentally dropping one ranking mode when sorting.
   Can rank all union candidates rather than impose a tiny final cap.
   Standard search.py rows: ["4x2",p,q,H,priority,ainvs].
3. Search.py checks full known fingerprint and parent duplicates, exact torsion,
   then ellrank(E,0); hunts effort=4 only if upper bound reaches target but too few
   points found. 45s timeout. Full verify_text before adding anything to board.
   A positive discovery writes board/first-...ec and research/discovery.json.
4. If none, possible extension toward full rank-7 height geometry:
   H >= log(p/2), so p <2*exp(12.948901), about 840k, suffices for all eligible
   models of this family. 600k covers much of it. Consider deeper prime bounds
   or more retained candidates to reduce heuristic misses.
5. On discovery stop unnecessary own searches, check entire board with parent
   research/frontier-board, require positive metrics.improvement. Full novelty,
   no duplicates/regression, total submission tree <=16MiB, all other caps.
6. Refresh frontier before upload, compare separately, only add inspected missing
   files or establish smaller held heights before replacing. No blind overwrites.
7. Public Yukon note must be >=5KiB and <=100KiB: objective, environment, sources,
   provenance, precise reproduction, actual failed experiments and results,
   measured gain, limitations. Strip secrets. Model GPT-6.1 Sol, harness Codex.
   Declare inherited/material dependencies explicitly: current frontier UUID above;
   family construction inspiration is carried in its inherited 3-r12 file.
   Do not claim inherited 104.611755 points as our gain.
   CLI submit --help and bundled yukon skill were already read.
   Claimed-score prefilter disabled.
8. Submit (user authorized), monitor official validation/promotion and verify
   actual user's leaderboard entry. Goal is not complete without positive entry.

## Space

At last check submission/ about 17MiB: clean generated intermediates before final
board verification / upload. Never remove board holders. Old generated files were
checksum-archived to /tmp/ec-rank-search-archive-20261025; manifest
research/archived-intermediates.json. Preserve source, winner, decisive logs and
provenance. Many older raw caches can likewise be archived safely after hash checks.
Research .ec proofs are ignored for score but count toward editable-tree cap.
No agents spawned; instructions do not authorize delegation.

## Approvals

Normal tools fail due WSL sandbox initialization; escalated WSL works.
Past automatic review rejected broad board replacement; separate comparison and
exclusive creation of only missing 4-9.ec resolved it. No unresolved approval
request. Earlier usage-limit failure was resolved when user resumed.
