#!/usr/bin/env bash
# Benchmark a submission.
#
#   1. Wipe the stale score.json before ANY fallible work, so a failed run
#      can never expose a result from an older invocation.
#   2. Resolve the venv and the pinned gp created by ./setup.sh; fail closed
#      if either is absent.
#   3. Run `python -E -s -m harness.verify` from the repo root (cwd on sys.path, so
#      the stage runs the checked-out ec_verify/ and harness/, the same tree the CLI
#      and the tests run; -E -s still ignore PYTHON* variables and the user site).
#      The harness executes NO competitor
#      code — the submission is a plain-text curve + points file — but it is
#      hostile input to our parser and to PARI/GP, so on Linux the verify
#      stage runs under bubblewrap (read-only filesystem, no network, no
#      capabilities, writable only in a throwaway scratch dir) and on macOS
#      under sandbox-exec. Without a sandbox we warn and run unconfined.
#   4. The harness writes score.json into the scratch dir (EC_SCORE_DIR); it
#      is copied to ./score.json only after success. Nonzero exit therefore
#      always means: no score file.
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${root}"

rm -f score.json

if [[ -x "${root}/.venv-bench/bin/python" ]]; then
  vpy="${root}/.venv-bench/bin/python"
else
  echo "!! venv missing; run ./setup.sh first" >&2
  exit 1
fi

PARI_VERSION="2.17.4"
gp_version_of() { "$1" --version-short 2>/dev/null | tr -d '[:space:]'; }
gp_bin=""
if [[ -x "${root}/tools/pari/bin/gp" ]]; then
  gp_bin="${root}/tools/pari/bin/gp"
elif command -v gp >/dev/null 2>&1; then
  gp_bin="$(command -v gp)"
fi
if [[ -z "${gp_bin}" || "$(gp_version_of "${gp_bin}")" != "${PARI_VERSION}" ]]; then
  echo "!! pinned gp ${PARI_VERSION} not found (got '${gp_bin:-none}' $(gp_version_of "${gp_bin:-/nonexistent}")); run ./setup.sh" >&2
  exit 1
fi

scratch="$(cd "$(mktemp -d -p "${EC_SCRATCH:-${TMPDIR:-/tmp}}")" && pwd -P)"
cleanup() { [[ -z "${scratch:-}" ]] || rm -rf "${scratch}" 2>/dev/null || true; }
trap cleanup EXIT
chmod 1777 "${scratch}"

# Parent board for the regression gate (docs/verifier-architecture.md §4). Yukon grades a single-parent commit on the
# previous frontier; benchmark.yml passes EC_PARENT_REF="${GITHUB_SHA}^" on submissions/* refs only
# (its EC_SUBMISSION_RUN flag; yukon/baseline/* and maintainer runs on main have no frontier),
# together with EC_REQUIRE_PARENT, which the harness honours: a run that reaches it without a
# parent board is rejected. Unset (local runs, baselines) => the gate is skipped; the score itself
# never depends on the parent. A failure before the harness speaks its interface too:
# `REJECTED: <CODE>: <detail>` on stdout, exit 1, no score.json.
parent_board=""
if [[ -n "${EC_PARENT_REF:-}" ]]; then
  # Attribute-free blob export (tools/export_parent_board.sh reads blobs directly), so a participant
  # .gitattributes cannot drop or rewrite a parent file. Bad ref or unreadable blob => rejection;
  # non-regular entries (symlinks, gitlinks, nested paths) are skipped with a `!!` note, as the
  # surface check already refuses them on submissions/* refs.
  parent_board="${scratch}/parent/submission/board"
  if ! bash "${root}/tools/export_parent_board.sh" "${root}" "${EC_PARENT_REF}" "${parent_board}"; then
    echo "REJECTED: BOARD_REGRESSION: the previous board could not be exported from EC_PARENT_REF '${EC_PARENT_REF}' (see the !! note above); the regression gate cannot run, which needs a maintainer"
    exit 1
  fi
  chmod -R a+rX "${scratch}/parent"
fi

run_verify=()
if command -v bwrap >/dev/null 2>&1; then
  bw=( bwrap )
  if [[ "$(id -u)" -ne 0 ]] && command -v sudo >/dev/null 2>&1 && sudo -n true >/dev/null 2>&1; then
    bw=( sudo -n bwrap )
  elif command -v setpriv >/dev/null 2>&1; then
    bw=( setpriv --no-new-privs bwrap )
  fi
  # Re-bind the repo at a world-traversable scratch-anchored path (a 0750
  # $HOME is not traversable once CAP_DAC_OVERRIDE is dropped).
  sandbox_repo="${scratch}/repo"
  mkdir -p "${sandbox_repo}"
  sandbox_vpy="${sandbox_repo}${vpy#"${root}"}"
  sandbox_gp="${gp_bin}"
  [[ "${gp_bin}" == "${root}/"* ]] && sandbox_gp="${sandbox_repo}${gp_bin#"${root}"}"
  run_verify=(
    "${bw[@]}"
      --ro-bind / / --dev /dev --ro-bind /proc /proc
      --bind "${scratch}" "${scratch}"
      --ro-bind "${root}" "${sandbox_repo}"
      --setenv TMPDIR "${scratch}"
      --setenv EC_SCORE_DIR "${scratch}"
      --setenv EC_GP_BIN "${sandbox_gp}"
      --setenv EC_PARENT_BOARD_DIR "${parent_board}"
      --setenv EC_REQUIRE_PARENT "${EC_REQUIRE_PARENT:-}"
      --chdir "${sandbox_repo}"
      --unshare-net --unshare-ipc --unshare-uts --unshare-cgroup
      --cap-drop ALL --new-session --die-with-parent
      -- "${sandbox_vpy}" -E -s -m harness.verify
  )
elif [[ "$(uname -s)" == "Darwin" ]] && command -v sandbox-exec >/dev/null 2>&1; then
  profile="(version 1)(allow default)(deny network*)(deny file-write*)(allow file-write* (subpath \"${scratch}\"))(allow file-write* (subpath \"/dev\"))"
  run_verify=(
    sandbox-exec -p "${profile}"
      /usr/bin/env TMPDIR="${scratch}" EC_SCORE_DIR="${scratch}" EC_GP_BIN="${gp_bin}" EC_PARENT_BOARD_DIR="${parent_board}" EC_REQUIRE_PARENT="${EC_REQUIRE_PARENT:-}"
      "${vpy}" -E -s -m harness.verify
  )
else
  echo "!! no sandbox available (bubblewrap/sandbox-exec); running verify UNCONFINED (dev fallback)" >&2
  run_verify=(
    env TMPDIR="${scratch}" EC_SCORE_DIR="${scratch}" EC_GP_BIN="${gp_bin}" EC_PARENT_BOARD_DIR="${parent_board}" EC_REQUIRE_PARENT="${EC_REQUIRE_PARENT:-}"
    "${vpy}" -E -s -m harness.verify
  )
fi

"${run_verify[@]}"

if [[ ! -s "${scratch}/score.json" ]]; then
  echo "!! harness exited 0 but produced no score.json" >&2
  exit 1
fi
cp "${scratch}/score.json" "${root}/score.json"
