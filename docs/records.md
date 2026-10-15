# The record store

How the known curves and the per-cell records were built and frozen. Written for
maintainers; participants need only `docs/submitting.md` §3. The files are frozen: a rebuild
changes the score of every board (§3). The scoring rule itself is `docs/verifier-architecture.md` §4.

## 1. Three files

Each layer is derived from the previous one by a tool in `tools/`. The verifier reads only
the second and third, from the read-only tree, so a commit on `main` that rebuilds them takes
effect for the next submission with no code change.

```text
tools/known_raw.json          raw mirror v2      every source record, with points and the certification result
ec_verify/known_curves.json   ledger v2          one compact entry per known curve: fingerprint, torsion, H, status
ec_verify/records.json        cell table         record of every (torsion, rank) cell; its height is the baseline (version 5)
```

The raw mirror and the ledger are written one record per line (`tools/fetch_known.dump_per_line`:
sorted top-level keys, then one compact JSON object per curve). They are still plain JSON for
`json.loads`, and a rebuild that adds a curve changes one line, so a pull request shows
exactly the curves it added. The cell table is small and uses `json.dumps(indent=1, sort_keys=True)`.

### 1.1 Raw mirror: `tools/known_raw.json`

Append-only, keyed by `id` (`<source>-<source_id>`). Entries that vanish upstream are kept.
Dujella rows are parsed from the cached pages `tools/cache/rk*.html` (rank records) and
`tools/cache/z*.html` (the per-torsion record pages and their "old records" pages) by
`tools/fetch_records.py:parse_curves`: every curve on every page with its independent points;
torsion points are not kept, every listed point is checked on the curve exactly, a misprinted point
is reported and left out, and a misprinted equation aborts the fetch unless the maintainer has
accounted for it in `fetch_known.DUJELLA_MISPRINTS` (one today: Jain's rank-31 curve on `z1.html`,
which is `icarm-3856`). The first curve of a rank page keeps the
id `dujella-rk20`; every other curve is `dujella-<page>-<8 hex of sha1(a-invariants)>`, stable
when a page is re-fetched with newer curves prepended. The LMFDB rows are the eight classic
anchors of `tools/anchor_curves.json` (no points) and `tools/lmfdb_curves.json`: the three
smallest curves by naive height of every (torsion, rank ≥ 1) in Cremona's ecdata tables
(conductor ≤ 500000, the LMFDB's curves over ℚ), with their generators and published torsion,
selected by `tools/fetch_lmfdb.py` from the `allgens` files (the LMFDB web API rate-limits bulk
queries behind a CAPTCHA). Ids are `lmfdb-<Cremona label>`.

```json
{"version": 2, "fetched_at": "2026-10-15T04:07:28+00:00", "curves": [
  {"id": "icarm-302", "source": "icarm", "source_id": "302",
   "ainvs": ["1", "1", "1", "-128...", "560..."],
   "points": [["-3437...", "-3100..."], ...],
   "rank_claimed": 31, "torsion_source": [], "naive_height_source": 468.277,
   "source_created_at": "2026-08-23 ...", "source_updated_at": "...",
   "first_seen": "2026-10-09T02:06:33+00:00",
   "verification": {"status": "certified", "certified_rank": 31, "code": null, "message": null,
                    "verified_at": "...", "pari_version": "2.17.4", "verifier_version": "1.0.0.dev0",
                    "source_updated_at": "..."}}
]}
```

`torsion_source`, `naive_height_source` (ICARM's `12·H`) and the `source_*_at` stamps are
`null` for sources that do not publish them. `verification` is written by
`tools/certify_known.py`; its `source_updated_at` records which version of the entry was
certified.

### 1.2 Ledger: `ec_verify/known_curves.json` (version 2)

One entry per raw entry, no points. Sorted by `rank_claimed` descending, then source (ICARM,
Dujella, LMFDB: a Dujella copy of an ICARM curve neither renames it in `KNOWN_CURVE` nor marks the
ICARM entry as the duplicate), then numeric id.

```json
{"version": 2, "frozen_at": "2026-10-15T04:07:28+00:00", "primes_upto": 1000, "curves": [
  {"id": "icarm-302", "source": "icarm", "source_id": "302",
   "ainvs": ["1", "1", "1", "-128...", "560..."],
   "fingerprint": [168 ints],
   "torsion": [], "torsion_key": "1", "naive_height": 39.023091,
   "rank_claimed": 31,
   "verification": {"status": "certified", "certified_rank": 31, "code": null,
                    "verified_at": "...", "pari_version": "2.17.4", "verifier_version": "1.0.0.dev0"},
   "duplicate_of_earlier": false,
   "source_created_at": "...", "first_seen": "..."}
]}
```

- `fingerprint` is `|a_p|` for the 168 primes below 1000, `−1` at bad primes. Every entry
  feeds the novelty gate, whatever its status.
- `verification.status` is `certified` (the published points certify; `certified_rank` may be
  below `rank_claimed`, and the curve then sits in the lower cell), `failed` (the points were
  rejected; `code` keeps the rejection code) or `unverified` (no points published). Only
  `certified` entries define records.
- `torsion` is computed by gp; `torsion_key` is `"1"`, `"n"` or `"nxm"`.
- `duplicate_of_earlier` marks an entry whose fingerprint already appeared earlier in the
  file (isogenous or re-uploaded curves); `known_curves.classes` groups entries by
  fingerprint and claimed rank.
- `frozen_at` is the raw mirror's `fetched_at`.

### 1.3 Cell table: `ec_verify/records.json`

```json
{"version": 5, "scoring_version": 6, "ledger_frozen_at": "2026-10-15T04:07:28+00:00",
 "built_at": "2026-10-15T...",
 "cells": {
   "1/6": {"torsion": [], "torsion_key": "1", "rank": 6,
           "record":   {"id": "icarm-56", "naive_height": 2.531337, "ainvs": ["0", "0", "1", "-277", "4566"],
                        "certified_rank": 6},
           "known_count": 1138, "class_count": 1108, "exact_count": 20}
 }}
```

A cell appears for every rank from 1 up to the highest certified rank of its torsion group;
every other cell is a breakthrough cell ("no known curve of that torsion reaches the rank").
A curve certified at rank R reaches every rank below R: list fewer of its points and the
verifier certifies the lower rank and places the curve in that cell. So the members of cell
(T, r) are all certified ledger curves of torsion T with certified rank at least r, and
`record` is the smallest of them (ties broken by the numerically smallest id). The record may
therefore be a higher-rank curve: `record.certified_rank` says which rank it was certified at,
and `exact_count` how many members sit at exactly r (0 for a rank no curve is certified at
exactly, which would otherwise have been a breakthrough cell worth more than 100). In the 2026-10-15 freeze
two cells, `2x2/13` and `3/12`, take their record from a curve one rank up. The record's height
is the cell's baseline (§3);
`known_count` counts ledger entries and `class_count` distinct fingerprints. `version`
increases by one on every rebuild; `scoring_version` must equal `SCORING_VERSION` in
`ec_verify/score.py` or the table does not validate. The verifier reads only `version` and
`cells.*.record`.

## 2. The pipeline (the freeze)

Run with the pinned gp from `./setup.sh`, in this order:

```bash
curl -fsSL https://web.math.pmf.unizg.hr/~duje/tors/z2.html -o tools/cache/z2.html   # ... and so on: refresh the cached Dujella pages by hand
python tools/fetch_lmfdb.py --ecdata /tmp/ecdata --download   # network: Cremona's 50 allgens files (180 MiB) -> tools/lmfdb_curves.json
python tools/fetch_known.py                 # network: ICARM database.json in one request; Dujella cache; LMFDB anchors + lmfdb_curves.json
python tools/certify_known.py --jobs 8      # certify new or updated entries (novelty off), written back into the raw mirror
python tools/build_known.py --jobs 8        # the ledger
python tools/build_records.py               # the cell table
python -m pytest tests -q                   # the acceptance check
```

The pipeline was run at the freeze of 2026-10-09 (table version 3: 4088 curves, 154 cells) and
again before launch on 2026-10-12 (version 4: 5541 curves, 156 cells; ICARM had grown by 126
curves, Dujella's per-torsion pages and the LMFDB selection were added, and the Z/2 cells 21 and
22 and 15 records moved), and once more on 2026-10-15 (version 5: 5586 curves, 156 cells). That run re-fetched ICARM (44 new
curves, all certified; one record moved, `4x2/7`) and re-cached Dujella's pages (one new curve, already on
ICARM); the LMFDB selection was not re-fetched, since Cremona's tables are static and `fetch_lmfdb.py` was
last run on 2026-10-12. Rebuilding changes the score of every board, so it is a maintainer
decision announced as a new `scoring_version` or `records_version`, never a scheduled job. All
the tools are deterministic and idempotent: run
twice on the same inputs, the second run changes nothing but timestamps and the table
`version`. `fetch_known.py` refuses malformed source records rather
than dropping them. `certify_known.py` takes `--force` to recertify everything and `--limit N`
for a trial run. `build_known.py` fails loudly if gp's torsion disagrees with a source's
published torsion. `build_records.py` fails loudly if a cell present in the previous table has
no certified curve any more (a certification that flipped, or an entry whose a-invariants
changed); that needs a maintainer, not a bot.

## 3. Frozen baselines

A cell's baseline is its record height and never changes. Curves published after the freeze
do not enter the table; they reach the board as ordinary submissions. Hence a board always
sums to the same number, and Yukon's stored
frontier score is always reachable.

## 4. No exemptions

The ledger is frozen, so a promoted curve can never become known later. The former
scoring-v5 `novelty_exempt` mark is gone. A parent-board file that matches the ledger therefore means the
ledger was rebuilt under the frontier: the run is rejected (`BOARD_REGRESSION`) until a
maintainer repairs the frontier board, rather than the cell being silently unheld.

## 5. Invariants the tests enforce

`python -m pytest tests -q` checks, over the shipped files:

- the ledger has every v2 field, a valid status, `certified_rank` exactly when certified, and
  168-entry fingerprints (`tests/test_known_curves_classes.py`, which also pins the curve and
  class counts: 5586 and 4239 in the 2026-10-15 freeze);
- every frozen curve matches its own fingerprint and has at least 128 good primes below 1000,
  the bound that keeps the candidate floor of 70 sound; the minimum is pinned at 139
  (`tests/test_novelty.py`);
- the cell table validates, its `scoring_version` is the code's, its `ledger_frozen_at` is the
  ledger's `frozen_at`, its cells are exactly the certified (torsion, rank) pairs of the
  ledger with the smallest curve that reaches each rank as record, no cell carries a `baseline`
  and the table no `h_cap`, all 15 torsion groups appear, and no torsion group has a rank gap below its top certified rank
  (`tests/test_records.py`);
- ICARM's published height is `12·H` to `10⁻⁵` on every raw entry that has one;
- on fixtures: the merge is idempotent and append-only, certification is incremental,
  fingerprints are reused, every rank up to a torsion group's top certified rank is a cell whose
  record is the smallest curve reaching it, non-integral published
  models are rejected rather than truncated, a re-fetched row whose points changed without a
  timestamp is re-certified, and no rank-0 cell is built
  (`tests/test_fetch_known.py`, `tests/test_tools_known.py`).

A rebuild that adds curves changes the pinned counts; update them in the same pull request.
If a new curve has fewer than 128 good
primes the novelty test fails, and `CANDIDATE_MIN_GOOD` in `ec_verify/fingerprint.py` must be
raised, which is a verifier change, not a data rebuild.

## 6. What a rebuild costs

Only new or changed entries cost anything. The ICARM fetch is one request; the merge keeps an
entry unchanged unless its `source_updated_at` or a-invariants moved, and a replaced entry
loses its `verification`, so `certify_known.py` runs only on new and updated ids. The initial
pass over the whole ICARM set takes about 15 minutes. `build_known.py` reuses
fingerprint and torsion from the previous ledger whenever id and a-invariants are unchanged;
`build_records.py` reads only the ledger.

## 7. How to add a source

Write rows in the raw v2 shape: `id` = `<source>-<source_id>`, a new `source` name,
`source_id` unique within it, `ainvs` as five integer strings, `points` as `[x, y]` string pairs
(`num/den` allowed; an empty list if the source publishes none, which makes the entry
`unverified`), `rank_claimed`, and `torsion_source`, `naive_height_source`,
`source_created_at`, `source_updated_at` (each `null` if unknown). Add a row function to
`tools/fetch_known.py` beside `dujella_rows` and `anchor_rows`, give the source a place in
`_SOURCE_ORDER` there (the build tools import it: a height tie between copies of one curve goes to
the source listed first, then to the numerically smallest id), and run the pipeline. If the source publishes torsion, `build_known.py` checks
it against gp. Nothing else changes: every ledger entry feeds the novelty gate, and those that
certify can define records.
