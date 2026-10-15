# Threat model

## Assets

- **A1 — board integrity:** no curve scores more than its points and its size prove, and
  no promoted record is lost.
- **A2 — determinism:** same bytes, same verdict, on any host.
- **A3 — runner confidentiality:** a submission cannot read or exfiltrate host files.
- **A4 — availability/cost:** a submission cannot burn unbounded time or memory.

## Attacker

The participant controls every byte under `submission/`: the number, names, types and
contents of the entries in `submission/board/`, and anything else beside it. They cannot
execute code on the runner, influence the checkout, the parent commit the workflow checks
out, the pinned PARI build, or the frozen record store (`ec_verify/known_curves.json`,
`ec_verify/records.json`).

## Controls

- **No participant code.** The only editable path holds data; the workflow admits only
  regular blobs there; the harness opens each board file `O_NOFOLLOW` and requires `S_ISREG`.
- **Board gate.** At most 512 `.ec` files (`TOO_MANY_CURVES`), names restricted to
  `[A-Za-z0-9._-]+.ec`, each file loaded strictly (regular file, 1 MiB cap, strict UTF-8)
  and run through every per-curve gate; fingerprints compared across files so a curve, its
  isogenies and twists count once (`DUPLICATE_CURVE`). The previous board for the regression
  gate is exported by `benchmark.sh` outside the sandbox, from the parent of the commit the
  workflow checked out, by `tools/export_parent_board.sh`: it reads each regular blob directly
  under `submission/board/` with `git cat-file blob` (no `git archive`, so no attribute
  processing: `export-ignore`, `export-subst` and filters do not apply), skips anything else
  with a `!!` note, and exits non-zero if any listed blob cannot be read (`set -euo pipefail`
  plus an explicit check on `git cat-file`), so a parent file is never silently missing.
  The workflow's surface check also rejects any dotfile under `submission/`
  (`.gitattributes`, `board/.gitattributes`) and any path with a character outside
  `[A-Za-z0-9._/-]`, so no odd name reaches the export.
- **Parser hardening.** ASCII only, `-?[0-9]+` tokens, digit/count/line/file caps, every
  failure a stable code. PARI scripts are assembled from validated integers, never from
  submission text.
- **Fail-closed acceptance.** Exact on-curve check; exact torsion stripping; independence
  only via the certified positive-definite pairing (§3 of the architecture). A tiny
  regulator is a rejection. `#CLAIM` is never trusted.
- **Pasting a public witness.** Every published high-rank curve (ICARM, Dujella) is in the
  frozen known set; the `|a_p|` fingerprint also rejects isogenies, twists and re-models
  (`KNOWN_CURVE`). Residual: a curve published elsewhere after the freeze, or in a source
  we did not ingest, is accepted by the gate. Submitting it is a rules violation
  (`RULES.md` §1–2), caught by the post-competition audit against ICARM, the tables, the
  literature and submission histories, with the reward revoked (§4).
- **Gaming the baseline.** The baseline is frozen in a non-editable file
  (`ec_verify/records.json`). `H` is computed from the minimal model's `c4` and `c6`, and the
  small-prime gate keeps models minimal; a model non-minimal at a large prime only has a
  larger `H`. The table never changes after the freeze, so a board always sums to the same
  number and Yukon's stored frontier score is always reachable (`docs/records.md` §3).
- **Sandbox.** Verify runs under bubblewrap (Linux) / `sandbox-exec` (macOS): no
  network, read-only tree, writable scratch only, capabilities dropped.
- **Budgets.** `parisizemax = 2 GiB`; per-stage timeouts and a 20-minute budget for the whole
  board (`RESOURCE_EXCEEDED`); 90-minute job timeout.
- **Score hygiene.** `score.json` wiped before any fallible work and copied out of
  scratch only on success; sorted keys; no timings.
- **Supply chain.** PARI tarball SHA-256 pinned (two independent downloads), setuptools
  pinned, SHA-pinned GitHub actions, `persist-credentials: false`, CI never runs a
  `submissions/**` head.

## Probed attacks

| attack | outcome |
|---|---|
| model rescaled by a prime above the 10⁶ minimality bound | undetected by the gate, but heights, torsion and the fingerprint are invariant and `H` grows: score can only fall |
| the rank-31 record, its re-model, a quadratic twist, an isogenous curve | `KNOWN_CURVE` |
| sextic twist of a known j = 0 curve | accepted — a different curve whose rank must be certified from scratch |
| `submission/.gitattributes` with `board/*.ec export-ignore` (empty or altered parent board) | blocked (workflow surface check rejects dotfiles under `submission/`) and ineffective (attribute-free parent export; an unreadable blob aborts it) |
| a submission graded on a ref the workflow does not recognise (no `EC_PARENT_REF`, gate skipped) | every `submissions/*` ref is gated against its parent commit, and `EC_REQUIRE_PARENT` makes the harness reject a run that reaches it without a parent board (`BOARD_REGRESSION`) |
| a file that makes gp error mid-script | gp runs with `recover = 0`, so the sentinel is never printed, and any non-warning `***` diagnostic on stderr is `GP_ERROR` naming the file |
| a valid 512-file frontier of worst-case files (31 points of ~1000 digits, or 20 of ~1800) to exhaust the next submission's budget | measured 0.8 s per file on an M-series laptop, about 9 minutes for 512 distinct files against the 1200 s budget; 64 dependent points of ~2000 digits reject in 11.5 s without a stack overflow |
| a frozen known curve with many bad primes below 200 (icarm-503, -319, -45, -374), pasted with its published points | `KNOWN_CURVE` (168 primes, fixed floor of 30 compared; self-match test over the whole freeze) |
| any of the 5586 frozen known curves, pasted with its published points | `KNOWN_CURVE` (the freeze is the whole ICARM database plus Dujella's pages and the LMFDB selection); the curves whose points certify passed `CURVE_NON_MINIMAL` in `certify_known`, so records use canonical models |
| quadratic twist by a smooth d of a known curve (icarm-51 by its 83 smallest good primes: 396-digit a4, minimal, 74 comparable primes) | `KNOWN_CURVE`: 74 ≥ 30; candidates with < 70 good primes below 1000 are `CURVE_TOO_MANY_BAD_PRIMES`, so ≥ 41 primes are always compared (architecture §3b) |
| generator plus its torsion translate P+T | `POINTS_DEPENDENT`, relation named |
| torsion point among generators | dropped exactly |
| 64 multiples of one generator, ~1500-digit coordinates | `POINTS_DEPENDENT`; certificate 7.6 s; hint path bounded to one 20 s attempt |
| rank-30 witness with every point replaced by 5·P (1181 digits) | certifies 30 in 1.3 s |
| 400-digit coefficients (1202-digit Δ) | every stage sub-second |
| `#CLAIM` with 5000 digits, directory / FIFO / symlink at the submission path | clean rejection |
| NUL, BOM, Unicode minus, `+6`, hex, exponent, signed denominator, inline comment | `PARSE_SYNTAX` |
| stage timeout | gp killed, no lingering process |

## Residual risk

- The certificate rests on `ellheight` accuracy (architecture §3). Mitigated by the
  two-precision tripwire; an independent recomputation (Sage or Magma) of the pairing matrix is
  required before announcing a rank-31-or-higher result (architecture §7) and is not performed
  for ordinary cell records.
- `metrics.curves[*].conductor` is wall-clock-bounded and may differ between machines; it
  is excluded from the determinism guarantee and never scored.
- ICARM rank claims are not trusted: only curves whose published points certify define
  records; a claimed record that fails stays visible as `failed`.
- The known set is what the freeze ingested (ICARM, Dujella's tables, a selection from the
  LMFDB). A curve published elsewhere, or after the freeze, passes the gate as new. Copying it
  is prohibited and reversed on audit, and the maintainers keep dated copies of the ICARM
  database for that comparison (`RULES.md` §4).
