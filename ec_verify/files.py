"""Board directory loading, shared by the harness and the local CLI (docs/verifier-architecture.md §5).

The board is data only: the loader refuses anything that is not a regular file named
`[A-Za-z0-9._-]+.ec` directly inside the directory, caps each file at MAX_FILE_BYTES and
decodes it as strict utf-8; check_submission_tree caps the whole editable tree at
MAX_BOARD_BYTES. Both entry points go through board.verify_board_dirs, so a board the laptop
accepts is a board the official run accepts.
"""

from __future__ import annotations

import os
import pathlib
import re
import stat

from .result import ErrorCode, VerifyError

MAX_FILE_BYTES = 1024 * 1024                    # per curve file (docs/verifier-architecture.md §5)
MAX_BOARD_BYTES = 16 * 1024 * 1024              # everything under submission/: benchmark.json maxSubmissionBytes, so the laptop and the platform agree
NAME = re.compile(r"[A-Za-z0-9_-][A-Za-z0-9._-]*\.ec")   # used with fullmatch; no leading dot (benchmark.yml rejects dotfiles)


def strict_load_text(path: pathlib.Path, max_bytes: int = MAX_FILE_BYTES) -> str:
    """Regular file only (no symlink, FIFO, device or directory), bounded read, strict utf-8.

    O_NOFOLLOW refuses a symlink at the final component; O_NONBLOCK makes opening a
    FIFO return instead of blocking for a writer; fstat/S_ISREG runs before any
    read so a directory or device is rejected rather than raising."""
    flags = (os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0) | getattr(os, "O_NONBLOCK", 0)
             | getattr(os, "O_BINARY", 0))
    try:
        fd = os.open(path, flags)
    except OSError as e:
        raise VerifyError(ErrorCode.SUBMISSION_NOT_REGULAR_FILE, f"cannot open {path.name}: {e}")
    try:
        if not stat.S_ISREG(os.fstat(fd).st_mode):
            raise VerifyError(ErrorCode.SUBMISSION_NOT_REGULAR_FILE, f"{path.name} is not a regular file")
        raw = b""
        while len(raw) <= max_bytes:
            chunk = os.read(fd, max_bytes + 1 - len(raw))
            if not chunk:
                break
            raw += chunk
    except OSError as e:
        raise VerifyError(ErrorCode.SUBMISSION_NOT_REGULAR_FILE, f"cannot read {path.name}: {e}")
    finally:
        os.close(fd)
    if len(raw) > max_bytes:
        raise VerifyError(ErrorCode.RESOURCE_EXCEEDED, f"{path.name} exceeds the {max_bytes}-byte cap")
    try:
        return raw.decode("utf-8", errors="strict")
    except UnicodeDecodeError as e:
        raise VerifyError(ErrorCode.PARSE_SYNTAX, f"{path.name} is not valid utf-8: {e}")


def load_board_dir(path: pathlib.Path, max_bytes: int = MAX_FILE_BYTES, on_error=None) -> dict:
    """{name: text} for every regular *.ec file directly inside `path`. A missing directory is an
    empty board. Non-.ec names and subdirectories are ignored; a .ec entry that is not a regular
    file, or whose name is outside [A-Za-z0-9._-], is a rejection — unless `on_error` is given,
    in which case it is called as on_error(name, code, message) and the entry is skipped (used
    for the parent board, whose unloadable entries must not reject the candidate)."""
    if path.is_symlink():                      # before exists(), which follows a dangling link to False
        raise VerifyError(ErrorCode.SUBMISSION_NOT_REGULAR_FILE, f"{path} is a symlink, not a directory")
    if not path.exists():
        return {}
    if not path.is_dir():
        raise VerifyError(ErrorCode.SUBMISSION_NOT_REGULAR_FILE, f"{path} is not a directory")
    out = {}
    for entry in sorted(os.scandir(path), key=lambda e: e.name):
        if entry.is_dir(follow_symlinks=False):
            continue
        if not entry.name.endswith(".ec"):
            continue
        try:
            if not NAME.fullmatch(entry.name):
                raise VerifyError(ErrorCode.PARSE_SYNTAX, "file name must match [A-Za-z0-9._-]+.ec with no leading dot")
            out[entry.name] = strict_load_text(pathlib.Path(entry.path), max_bytes)
        except VerifyError as e:
            if on_error is None:
                raise VerifyError(e.code, f"board/{entry.name}: {e.message}")      # the documented detail prefix
            on_error(entry.name, e.code.value, e.message)
    return out


def load_board_dir_lenient(path: pathlib.Path, max_bytes: int = MAX_FILE_BYTES) -> tuple:
    """load_board_dir for the parent board: entries that cannot be loaded are skipped and
    returned as {file, code, message} records (metrics.ignored_parent_files) rather than
    rejecting the candidate. Returns (files, ignored)."""
    ignored: list = []
    files = load_board_dir(path, max_bytes, on_error=lambda name, code, message:
                           ignored.append({"file": name, "code": code, "message": message}))
    return files, ignored


def submission_tree_bytes(path: pathlib.Path) -> int:
    """Bytes of every regular file under `path`, recursively, symlinks not followed and not counted:
    what the platform's maxSubmissionBytes measures over the editable tree. A missing path is 0, and
    so is a file that vanished between listing and stat (a search loop rewriting scratch under
    submission/). Anything else that cannot be read (an unreadable directory, EACCES on a file) is a
    rejection: an undercount here would be a submission the platform refuses."""
    def unmeasurable(e: OSError) -> None:
        if isinstance(e, FileNotFoundError):            # the root is missing, or an entry vanished mid-walk
            return
        raise VerifyError(ErrorCode.RESOURCE_EXCEEDED,
                          f"{e.filename}: cannot be read ({e.strerror}); the submission tree could not be measured")

    total = 0
    for dirpath, dirnames, filenames in os.walk(path, onerror=unmeasurable, followlinks=False):
        for fn in filenames:
            try:
                st = os.lstat(os.path.join(dirpath, fn))
            except OSError as e:
                unmeasurable(e)
                continue
            if stat.S_ISREG(st.st_mode):
                total += st.st_size
    return total


def check_submission_tree(path: pathlib.Path, max_bytes: int | None = None) -> int:
    """RESOURCE_EXCEEDED when the editable tree is over the platform's cap, so the laptop refuses a
    submission the platform would refuse (benchmark.json maxSubmissionBytes covers all of submission/,
    not only board/*.ec). Called once, from board.verify_board_dirs. Returns the byte count."""
    if max_bytes is None:
        max_bytes = MAX_BOARD_BYTES
    total = submission_tree_bytes(path)
    if total > max_bytes:
        raise VerifyError(ErrorCode.RESOURCE_EXCEEDED,
                          f"{path.name}/ totals {total} bytes; at most {max_bytes} (the platform's submission cap "
                          "covers everything under submission/, not only the board)")
    return total
