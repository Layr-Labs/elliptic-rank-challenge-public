# Hermes EC Rank helper

This directory is an independent, minimal helper for local EC Rank work. It
does not alter the verifier, the frozen record table, or the promoted board.

`random_rank_search.py` performs a small deterministic/random search over short
Weierstrass models and uses the pinned PARI/GP binary from `tools/pari/bin/gp`
to screen curves by torsion and Mordell-Weil rank. It is intentionally a
screening tool, not a certificate generator: a curve is only eligible for a
board after the official `ec_verify` harness has checked minimality, novelty,
exact points, torsion, and height-pairing independence.

Usage:

```sh
python3 submission/search/hermes-search/random_rank_search.py \
  --count 100 --rank 1 --seed 1
```

The script prints candidate `a1 a2 a3 a4 a6` tuples. It does not write into
`submission/board` and does not claim a record.
