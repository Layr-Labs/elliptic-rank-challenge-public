#!/usr/bin/env bash
# Set up the benchmark environment. This is the ONLY step that may use the
# network. After it succeeds, ./benchmark.sh runs fully offline.
#
#   1. Python >= 3.11 venv at ./.venv-bench: an interpreter only. The verifier
#      is stdlib-only and benchmark.sh runs the checked-out tree from the repo
#      root (`python -E -s -m harness.verify`), so nothing is installed into
#      the venv: one copy of the code, no snapshot to go stale. benchmark.sh
#      re-derives the venv from that absolute path — exported variables do NOT
#      persist between Yukon's setupCommand and benchmarkCommand.
#   2. PARI/GP pinned to ${PARI_VERSION}: a gp on PATH that reports exactly
#      that version is accepted (dev machines: Homebrew); otherwise the
#      official tarball is downloaded, sha256-checked, and built into
#      tools/pari. benchmark.sh resolves the binary the same way
#      (tools/pari/bin/gp first, then PATH) and passes it as EC_GP_BIN.
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MIN_MINOR=11
PARI_VERSION="2.17.4"
PARI_URL="https://pari.math.u-bordeaux.fr/pub/pari/unix/pari-${PARI_VERSION}.tar.gz"
# sha256 of the official tarball (recorded from two independent downloads on 2026-09-09).
PARI_SHA256="02651d99c391007d384b3fadbc20abc6916b77036f9e496c99e9ce8688ca4b53"

find_python() {
  local candidate
  for candidate in "${PYTHON:-}" python3.13 python3.12 python3.11 python3 python; do
    [[ -n "${candidate}" ]] || continue
    command -v "${candidate}" >/dev/null 2>&1 || continue
    if "${candidate}" -c "import sys; raise SystemExit(0 if sys.version_info >= (3, ${MIN_MINOR}) else 1)" 2>/dev/null; then
      command -v "${candidate}"
      return 0
    fi
  done
  return 1
}

python_bin="$(find_python)" || {
  echo "setup.sh: no Python >= 3.${MIN_MINOR} found (set PYTHON=/path/to/python3)" >&2
  exit 1
}
echo "python:  ${python_bin} ($("${python_bin}" -c 'import sys; print(sys.version.split()[0])'))"

venv="${root}/.venv-bench"
if [[ ! -x "${venv}/bin/python" ]]; then
  "${python_bin}" -m venv "${venv}"
fi
vpy="${venv}/bin/python"

# --- PARI/GP -----------------------------------------------------------------
gp_version_of() { "$1" --version-short 2>/dev/null | tr -d '[:space:]'; }

gp_bin=""
if [[ -x "${root}/tools/pari/bin/gp" && "$(gp_version_of "${root}/tools/pari/bin/gp")" == "${PARI_VERSION}" ]]; then
  gp_bin="${root}/tools/pari/bin/gp"
elif command -v gp >/dev/null 2>&1 && [[ "$(gp_version_of "$(command -v gp)")" == "${PARI_VERSION}" ]]; then
  gp_bin="$(command -v gp)"
fi

if [[ -z "${gp_bin}" ]]; then
  echo "pari:    building pari-${PARI_VERSION} from source (no pinned gp found)"
  src="${root}/tools/pari-src"
  mkdir -p "${src}"
  tarball="${src}/pari-${PARI_VERSION}.tar.gz"
  if [[ ! -s "${tarball}" ]]; then
    curl -fsSL --retry 3 -o "${tarball}" "${PARI_URL}"
  fi
  if command -v sha256sum >/dev/null 2>&1; then
    echo "${PARI_SHA256}  ${tarball}" | sha256sum -c - >/dev/null
  else
    echo "${PARI_SHA256}  ${tarball}" | shasum -a 256 -c - >/dev/null   # macOS without coreutils
  fi
  rm -rf "${src}/pari-${PARI_VERSION}"
  tar -xzf "${tarball}" -C "${src}"
  (
    cd "${src}/pari-${PARI_VERSION}"
    # --static: gp links libpari statically, so the sandboxed harness (minimal environment,
    # repo re-bound at another path) needs no loader configuration.
    ./Configure --static --prefix="${root}/tools/pari" --without-readline >/dev/null
    make -j"$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 2)" gp >/dev/null
    # install-bin-sta installs gp-sta as gp; the docs targets need TeX and are skipped.
    make install-bin-sta install-cfg install-data >/dev/null 2>&1 || make install-bin-sta >/dev/null
  )
  gp_bin="${root}/tools/pari/bin/gp"
  [[ "$(gp_version_of "${gp_bin}")" == "${PARI_VERSION}" ]] || { echo "!! built gp reports $(gp_version_of "${gp_bin}")" >&2; exit 1; }
fi
echo "pari:    ${gp_bin} ($(gp_version_of "${gp_bin}"))"

# Optional: bubblewrap gives benchmark.sh its no-network sandbox on Linux.
if [[ "$(uname -s)" == "Linux" ]] && ! command -v bwrap >/dev/null 2>&1; then
  SUDO=""
  if [[ ${EUID:-$(id -u)} -ne 0 ]] && command -v sudo >/dev/null 2>&1 && sudo -n true >/dev/null 2>&1; then
    SUDO="sudo -n"
  fi
  if command -v apt-get >/dev/null 2>&1 && [[ -n "${SUDO}" || ${EUID:-$(id -u)} -eq 0 ]]; then
    ${SUDO} apt-get install -y --no-install-recommends bubblewrap || true
  fi
fi

echo "venv:    ${vpy}"
echo "harness: $(cd "${root}" && "${vpy}" -E -s -c 'import ec_verify; print(ec_verify.__version__)')"
