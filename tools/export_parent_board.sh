#!/usr/bin/env bash
# Export <ref>:submission/board into <dest> for the board regression gate (docs/verifier-architecture.md §4).
#
#   usage: tools/export_parent_board.sh <repo-root> <ref> <dest-dir>
#
# Attribute-free by construction: blobs are read with `git cat-file blob`, so no
# .gitattributes (export-ignore, export-subst, filters) inside the participant-editable
# tree can drop or rewrite a parent file. Only regular blobs (100644/100755) directly
# under submission/board/ are written; anything else (nested paths, symlinks 120000,
# gitlinks 160000) is skipped with a `!!` note; a blob that cannot be read exits 1. A ref
# without submission/board yields an empty <dest>; a ref that is not a commit exits 1.
set -euo pipefail

export_parent_board() {
  local root="$1" ref="$2" dest="$3"
  if ! git -C "${root}" rev-parse --verify --quiet "${ref}^{commit}" >/dev/null; then
    echo "!! EC_PARENT_REF '${ref}' is not a commit in this checkout (fetch-depth >= 2 required)" >&2
    return 1
  fi
  mkdir -p "${dest}"
  # Capture the listing first so a git ls-tree failure is checked (a process substitution
  # would swallow it); NUL-separated entries cannot be held in a shell variable. The file
  # sits beside <dest> (writable by construction), never inside it, and is removed after.
  local listing
  listing="$(mktemp "${dest%/}.ls-tree.XXXXXX")"
  if ! git -C "${root}" ls-tree -r -z "${ref}" -- submission/board > "${listing}"; then
    rm -f "${listing}"
    echo "!! parent board: git ls-tree failed for ${ref}" >&2
    return 1
  fi
  # Entries are NUL-separated: a blob name may contain a newline.
  local prefix="submission/board/" rec meta path mode type sha name
  while IFS= read -r -d '' rec; do
    meta=${rec%%$'\t'*}
    path=${rec#*$'\t'}
    read -r mode type sha <<< "${meta}"
    name=${path#"${prefix}"}
    if [[ "${path}" != "${prefix}"* || -z "${name}" || "${name}" == */* ]]; then
      echo "!! parent board: skipping nested or out-of-board entry ${path} (mode ${mode})" >&2
      continue
    fi
    if [[ "${type}" != "blob" || ( "${mode}" != "100644" && "${mode}" != "100755" ) ]]; then
      echo "!! parent board: skipping non-regular entry ${path} (mode ${mode})" >&2
      continue
    fi
    if ! git -C "${root}" cat-file blob "${sha}" > "${dest}/${name}"; then
      rm -f "${listing}"
      echo "!! parent board: cannot read blob ${sha} for ${path}" >&2
      return 1
    fi
  done < "${listing}"
  rm -f "${listing}"
  chmod -R a+rX "${dest}"
}

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
  if (( $# != 3 )); then
    echo "usage: $0 <repo-root> <ref> <dest-dir>" >&2
    exit 2
  fi
  export_parent_board "$@"
fi
