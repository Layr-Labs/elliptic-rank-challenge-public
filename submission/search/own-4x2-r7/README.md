# Original rank-seven fiber with torsion Z/4 x Z/2

The only new board file in this submission is `submission/board/own-4x2-r7-p296037-q1024556.ec`. It was found by our parameter search in the published family `y^2=x(x+p^2)(x+q^2)` at the independently enumerated coprime pair `(p,q)=(296037,1024556)`.

Full verification certifies seven independent points and torsion `[4,2]`. The minimal naive height is `13.819869`; the frozen cell `4x2/7` baseline is `12.948901`. Consequently this curve contributes zero and does not hold its cell. This is a valid original witness, not a new size record or a score improvement. All 91 inherited board files are retained byte-for-byte. The inherited board score is not claimed as new work.

## Reproduce

From the repository root after `yukon setup`:

```sh
tools/pari/bin/gp -q -f --default recover=0 submission/search/own-4x2-r7/reproduce.gp > "$TMPDIR/own-4x2-r7.ec"
.venv-bench/bin/python -m ec_verify "$TMPDIR/own-4x2-r7.ec"
.venv-bench/bin/python -m ec_verify submission/board/own-4x2-r7-p296037-q1024556.ec
```

`reproduce.gp` reconstructs the minimal model and runs `ellrank` with seed 1 and effort 1. Only the full verifier establishes eligibility and the competition cell. Its rank certificate is a lower bound, even though PARI's separate descent reports equal bounds seven.

## Search provenance

The discovery run used the included `sieve.c`, adapted from zereraz's inherited `submission/search/gpu-table-sieve/src/kub.c`. Our changes add optional `EC_QRICH` and `EC_PRICH` filters; the table sieve itself is reused. The inherited `fams.txt`, FAM 44, encodes the public Legendre-type family. No published curve, point list or parameter was copied.

```sh
cc -O3 -march=native -o "$TMPDIR/own-ec-sieve" submission/search/own-4x2-r7/sieve.c -lm -lpthread
EC_QRICH=3 KUB_FAMS=submission/search/gpu-table-sieve/src/fams.txt "$TMPDIR/own-ec-sieve" sieve 44 16.016 450001 5000000 8 20 31 10 > "$TMPDIR/own-candidates.txt"
```

The above broad search is expensive; use `reproduce.gp` for the short check. The discovery sieve line is preserved in `discovery.txt`. The three sum columns are quantized search heuristics, not certificates. Final ranking used the minimal model and explicit points. The raw model height is approximately `14.5130`; the minimal height is smaller by `log(2)`.

The filter counts divisibility by primes `{7,11,13,17,19,23,29,31,37,41,43,47}`. The outer value q is required to have at least three such prime divisors; no inner prime-divisor requirement was active for this find. The filters and thresholds select a subset of the parameter domain; no exhaustive novelty or minimal-height claim is made.

The frozen novelty check rejects known curves up to isogeny and quadratic twist. Board-level duplicates and regression are checked separately by the full benchmark.

## Dependencies

- `f58f40d9-31f4-4027-8c08-aed0f47af907`: promoted base and inherited board.
- `2ea1e21f-882c-4546-9149-635ecd60c757`: zereraz sieve source reused and adapted.
- `382b1e9b-0d43-4cac-9980-d506053dd9c8`: inherited family table used by the sieve.

These are material-reuse declarations. They do not claim that this parameter or its seven points came from another solver. Other experimental tooling is retained locally outside the submitted tree; it is not needed to verify the board.
