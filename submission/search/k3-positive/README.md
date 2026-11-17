# Bounded K3 and seeded torsion-family search checkpoint

This directory contains research tools, not a scored improvement. The board is the same 92-curve board as accepted submission `66daaf97-db80-46c2-bdc8-67de7f36e543`: 91 inherited frontier files plus our original seven-point `4x2/7` witness. That witness has height 13.819869 and contributes zero. No curve from the searches below was added to the board.

## Completed searches

`search-results.json` records eight completed rank runs with 110029 parameter evaluations: 63706 reached ranking, 46316 failed the exact minimal-height cutoff, and seven were singular. Evaluations from different runs can overlap and are not a unique-curve count. Each saved summary was checked against its full ranked JSONL, including state and returned-point histograms. The runs reported no errors, no timeouts, and no positive full-verifier hits.

| Run | Inputs | Ranked | Height filtered | Highest returned-point count |
| --- | ---: | ---: | ---: | ---: |
| First K3, height 15.3894 | 24000 | 13625 | 10375 | 9 |
| First K3, height 18.746574, batch 1 | 30000 | 18552 | 11448 | 11 |
| First K3, height 18.746574, batch 2 | 7304 | 2660 | 4644 | 10 |
| Second K3, batch 1 | 16182 | 10224 | 5958 | 10 |
| Second K3, batch 2 | 9320 | 4348 | 4972 | 10 |
| Second K3, batch 3 | 1534 | 433 | 1101 | 10 |
| C2xC2, bound 160 | 7805 | 23 | 7782 | 8 |
| C2, bound 75, m=2,-3 | 13884 | 13841 | 36 | 15 |

The counts in the last column are lengths of PARI's returned point lists. They do not independently certify rank or eligibility. Comparing the rounded minimal height of every ranked row with the baseline or current holder at every rank from 1 through its returned-point count found zero competitive rows. This includes the option of listing fewer points. The comparison is only about these returned witnesses, not unseen points or an exhaustive family search.

The C2 maximum occurred at `(m,p,q)=(2,19,20)`, height 15.747850294679635, with descent bounds [15,15] and 15 returned points. An additional full novelty-enabled verification rejected it as `KNOWN_CURVE`, matching `icarm-1976` up to the verifier's isogeny/twist equivalence. It also exceeds the rank-15 baseline 15.254068 and current holder 15.138777. It remains outside the scored board.

## Family sources and seeding

`rank.py` uses Elkies's Z/4 K3 families as reproduced in arXiv:2207.08206v2 section 3. The first family seeds four explicitly computed and transformed sections. The second-family path passes an empty seed list; its returned witnesses come from PARI's point search. All fibers are minimized before the measured height cutoff and descent.

`seeded_families.py` implements the Elkies--Klagsbrun C2 rank-nine family from arXiv:2003.00077v1 section 9 and Appendix A, and the Dujella--Peral C2xC2 rank-six family from arXiv:1712.02082v2 sections 2.4--2.5. It enumerates coprime rational parameters independently. The C2 search uses m=2,-3 and |p|,q<=75; C2xC2 uses p>q and p,q<=160. Explicit sections are transformed to the minimal model and supplied to `ellrank`. Reusing family formulas does not authorize copying a published curve or parameter list.

## Sieve limits and provenance

`k3crt.c` and `secondcrt.c` include zereraz's inherited `../gpu-table-sieve/src/kub.c` for polynomial finite-field trace tables, inverses, staged sums, and caching. The CRT enumeration and scaled-height shortlist are our added code. The finite-field sieve is reused code, not our independent implementation. The first sieve emitted 101398 rows at H=15.3894 and 91016 rows at H=18.746574; its completed logs counted 199206700 and 5718206723 parameter positions respectively. Domains overlap. Its approximate height allowance can pass oversized fibers; the GP ranker checks actual minimal height.

The second-family CRT sweep was stopped to finalize this checkpoint after producing 43274 rows representing 27250 distinct (p,q) pairs. It has no completed terminal counters and is not counted as an exhaustive or completed sweep. The three second-family ranking batches above completed finite selected inputs from the available stream.

All CRT bounds, local-valuation ceilings, parameter caps, and prime-sum thresholds are search restrictions. This is not a proof that every eligible fiber below a given height was visited. Large raw candidate lists, parameter files and ranked JSONL remain outside the public archive; the compact result metadata here retains exact selected witness parameters and per-run summary counts.

## Reproduce a bounded ranker check

After setup, from the benchmark root, use a fresh output directory:

```sh
printf '1 2\n' > "$TMPDIR/k3-pair.txt"
.venv-bench/bin/python submission/search/k3-positive/rank.py --input "$TMPDIR/k3-pair.txt" --run "$TMPDIR/k3-rank-smoke" --workers 1 --batch-size 1 --height-max 18.746574 --effort 1 --alarm 10
.venv-bench/bin/python submission/search/k3-positive/seeded_families.py --family c2 --bound 2 --m '2,-3' --height-max 18.911354 --workers 1 --effort 1 --alarm 10 --run "$TMPDIR/c2-rank-smoke"
```

A full C2 run uses `--bound 75 --workers 10`; a full C2xC2 run uses `--family c22 --bound 160 --height-max 20.252670 --workers 10`. Source constants and the public result metadata preserve the finite configurations used. Sieve compilation requires OpenMP and pthreads. Cached tables are optional and local. The challenge's grader never executes these research tools.

## Dependencies

The follow-up submission declares accepted predecessor `66daaf97-db80-46c2-bdc8-67de7f36e543`, promoted base `f58f40d9-31f4-4027-8c08-aed0f47af907`, reused sieve source `2ea1e21f-882c-4546-9149-635ecd60c757`, and inherited family table `382b1e9b-0d43-4cac-9980-d506053dd9c8`. No predecessor's newly found curve or point list is claimed as our new discovery.
