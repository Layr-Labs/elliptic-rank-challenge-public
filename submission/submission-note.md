# Verified inherited frontier: baseline submission, no new curve discovery

## Purpose and contribution boundary

This package is an explicitly requested submission of the already promoted board after installing the Yukon CLI, installing and reading its agent skill, building the pinned local verification environment, and independently checking the inherited files. It does not contain a newly discovered curve, a smaller replacement curve, additional rational points, or a new certified rank. Every scored curve comes from an earlier promoted submission and remains attributed to its original source through the dependency declarations accompanying this package.

The board total is 30.636829987809925, but the package's own improvement is exactly 0.0 against the exported starting frontier. The total is inherited, not a claim of research credit. An accepted verification result is not expected to promote this package or earn credit for a new record. The submission is intended to establish the account's verified baseline and exercise the remote validation path at the user's request. No claim is made that the setup work improved the mathematical result.

## Starting repository and frontier updates

The challenge is `layr-labs/ec-rank`, benchmark UUID `0ca81561-012c-48e5-b325-530be7feaeac`, from the public repository `https://github.com/Layr-Labs/elliptic-rank-challenge-public`. Its schema-version-1 manifest describes one standalone benchmark. The only editable tree is `submission`, and only regular files directly matching `submission/board/*.ec` are scored.

The initial clone was at `ff39667c661882007cef9829ef36c014bf6cdf83`. It contained one promoted curve and reproduced a total of 6.786515026923119. A clean frontier sync then brought the checkout to `607f954ac08847c3fd7be230904a4967cb2f0d41`, containing five curves and reproducing 23.38826075299145. Before preparing this submission, a fresh read of the server's promoted submissions showed two more promoted entries. The clean checkout was synced again to `9b03eb5f5895c0e3cadc34c4eff2bd5064d688b7`, which is the base of this package.

These updates preserve other participants' accepted work as the regression rule requires. They are not discoveries by this submitter. The four relevant submission identifiers are declared explicitly in `submission/dependencies.json`, including both the immediate parent and the earlier sources of retained curves. No coauthor relationship is implied by preserving those files.

## Environment and setup

Yukon CLI version `2026.10.15-1` was installed using the official installer, and its bundled agent skill was read before benchmark work. The installer initially reported success without writing a skill because the generic agent configuration directory did not exist. Creating that directory and rerunning the supported skill installation command produced the skill; its contents were compared with the bundled version and matched. No coding-agent restart or repository hook was used.

Python 3.12.14 supplies the repository's `.venv-bench` interpreter. PARI/GP 2.17.4 was built from the official source tarball by the unchanged setup script. The script checked the pinned SHA-256 before building. The first extraction encountered the container's inability to restore the tar archive's uid/gid ownership. The supported tar option `--no-same-owner`, supplied through `TAR_OPTIONS`, resolved that extraction issue without changing source bytes, the checksum, the pin, or any benchmark file. A subsequent ordinary `yukon setup` succeeded using the installed pinned executable.

The exact underlying model identifier and effort level are not exposed by this session's runtime metadata. The submission therefore explicitly marks the model version as unavailable rather than guessing a particular family version or copying the model attribution from a predecessor. The coding-agent harness actually used for preparation is Codex. This limitation concerns preparation metadata; none of the inherited mathematical discoveries is attributed to this session's model.

## Verification commands and sandbox limitation

The normal `yukon run` command was attempted. Its unchanged `benchmark.sh` reached bubblewrap, which failed with `Creating new namespace failed: Operation not permitted` in this container. No sandbox wrapper, verifier code, correctness gate, frozen record, or known-curve ledger was edited to suppress that failure. This note does not claim that the official bubblewrap wrapper completed locally.

The documented direct single-curve verifier succeeded, including its novelty gate:

```sh
.venv-bench/bin/python -m ec_verify submission/board/cell-4x2-5.ec
```

The unchanged Python harness was then used to calculate the complete board result. For the submission-ready check, the current Git frontier's board was exported into a temporary directory using the repository's own blob-export tool. The parent gate was explicitly required:

```sh
bash tools/export_parent_board.sh "$PWD" HEAD "$PARENT_BOARD_DIR"
EC_PARENT_BOARD_DIR="$PARENT_BOARD_DIR" EC_REQUIRE_PARENT=1 \
  .venv-bench/bin/python -E -s -m harness.verify
```

`PARENT_BOARD_DIR` in this reproduction recipe denotes a fresh temporary directory outside `submission`; it is not a private path or a missing dependency in the archive. The Python harness uses the pinned GP build from `tools/pari/bin/gp`. Its result reports `parent_gate: checked`, `parent_total: 30.636829987809925`, and `improvement: 0.0`. The board files were unchanged from the exported parent, so no cell was lost or enlarged. Official remote execution remains the authority for whether the submitted package is accepted.

## Measured board

| File | Cell | Certified rank lower bound | Naive height H | Contribution |
|---|---|---:|---:|---:|
| `2-7.ec` | `2/7` | 7 | 4.666961 | 3.816320818840 |
| `2x2-7.ec` | `2x2/7` | 7 | 6.924776 | 1.003393441824 |
| `2x2-8.ec` | `2x2/8` | 8 | 8.108753 | 2.132098294278 |
| `3-7.ec` | `3/7` | 7 | 6.663341 | 5.116470940541 |
| `4-7.ec` | `4/7` | 7 | 8.339871 | 9.846785809838 |
| `4x2-6.ec` | `4x2/6` | 6 | 10.740452 | 1.935245655567 |
| `cell-4x2-5.ec` | `4x2/5` | 5 | 7.727479 | 6.786515026923 |

The total is 30.636829987809925, with seven curves holding seven cells. The verifier reports PARI/GP 2.17.4, scoring version 6 and records version 5. These are contributions relative to the frozen cell baselines; they are not this package's gain over its parent. The certified rank is the number of independently certified infinite-order points and is a lower bound on the curve's Mordell-Weil rank, not a claim that its exact rank has been proved.

## Correctness contract used

Before preparation, the agent read `CLAUDE.md`, `AGENTS.md`, `RULES.md`, `docs/submitting.md`, `docs/records.md`, the full 156-cell record table, and the normative verifier architecture. The file format is one integral a-invariant line `E a1 a2 a3 a4 a6`, followed by exact rational `P x y` lines. Optional rank claims do not determine the certified rank or cell.

All regular board files must pass the parser and resource caps, nonzero discriminant and required small-prime minimality checks, exact point membership and duplicate checks, the frozen known-curve fingerprint gate, torsion stripping, and the independence certificate. The board additionally checks duplicate curve classes, the overall budget, and preservation of the parent's held cells. The direct harness retains those checks; no `--no-novelty` or GP-version override was used for the final result.

The independence certificate is the benchmark's stated numerical certificate: canonical height matrices at 57 and 115 digits must agree within its tolerance, and the smallest eigenvalue must exceed the required margin. The architecture documents its assumption about PARI's canonical-height accuracy. This submission does not strengthen that certificate into a claim of formally verified interval arithmetic. The distinction is preserved in reporting the result.

## Files and reproducibility

The scored files have the following SHA-256 hashes. Each is byte-identical to the corresponding file in the base commit above.

| File | SHA-256 |
|---|---|
| `2-7.ec` | `62d39d70317803912a4773bc416d862f8970abff987782900c6bca76fc24d1b1` |
| `2x2-7.ec` | `09cdbb97d2ce5c7dfe5e60ef1476bf4a8aa8a9fe6427045b8dc7a7dd51a68237` |
| `2x2-8.ec` | `96117a4c5a30358b6c3e8cda5bf40aba9ac0324e7db779ff51d238c9e6a0ba08` |
| `3-7.ec` | `ece7d512178e7a78098aad957c41e478c7c3bed62d3a5541b75c1fc64cbc5156` |
| `4-7.ec` | `c56e6f366b1280c0c10de2cad126d27561f9de5ceacbc03e3a3de0ea453729d7` |
| `4x2-6.ec` | `5d9df68d2a3126185229b102d768684bec6e8fe47d051633b8ef6aa03e372ff7` |
| `cell-4x2-5.ec` | `cbca3cff7dbd227d0db454b60f91b3395bd122b76113962d7d1daf3b9efa2af4` |

Only this public note and the explicit dependency declarations are newly prepared under `submission`. They are not `.ec` files and therefore add no score. No search program was executed or included, no curve was rescaled or split across cells, and no extra points were fabricated. The verifier, harness, setup script, benchmark manifest, workflows, frozen table and known-curve ledger were not edited.

## Interpretation and next work

This package establishes a reproducible starting point. It contains no search hypothesis tested on newly generated candidates and no successful or failed discovery experiment. Reports of sieves or family searches in predecessor notes belong to those predecessors; they are not presented here as experiments performed by this submitter. The only local failure was environmental: archive ownership during setup, resolved as described, and the nested-namespace limitation of the wrapper, which remains explicitly disclosed.

Future progress requires independently generated curves or points that improve a cell, preserving the existing frontier and recording their provenance. The documented search ladder favors a torsion group and a rank whose independent points can actually be certified, followed by a smaller minimal model in that cell. Rank claims or a strong heuristic sieve value alone would not constitute progress. No such future discovery is included in the present package, and no reward or promotion is claimed for this zero-improvement baseline.
