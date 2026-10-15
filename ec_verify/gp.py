"""The one place PARI/GP is invoked.

A script is assembled from *validated integers only* (never from submission
text), written into the scratch directory, and run as `gp -q -f --default=parisize=…
--default=parisizemax=… --default=recover=0 script` under a per-stage timeout, clipped to the run's deadline when one is set. Results come back as `KEY=VALUE` lines; a
sentinel line proves the script ran to completion, and any non-warning `***`
diagnostic on stderr is an error as well (with `recover=0` gp aborts on any error, so
the sentinel is missing; the stderr check names the error instead of an anonymous failure). GP prints nothing for a
one-line script without a trailing newline, so scripts are always
multi-line and end in `quit;\\n`.
"""

from __future__ import annotations

import contextlib
import contextvars
import os
import pathlib
import re
import shutil
import subprocess
import tempfile
import time

PINNED_VERSION = "2.17.4"
PARISIZE = 2 ** 28               # 256 MiB initial stack; gp grows it (with a stderr warning) up to PARISIZEMAX on large inputs
PARISIZEMAX = 2 ** 31            # 2 GiB stack ceiling (big-integer heights on 100-digit curves)
DEFAULT_PRECISION = 38           # decimal digits of realprecision unless a stage asks for more
SENTINEL = "EC_DONE"

_INT_RE = re.compile(r"^-?[0-9]+$")
_GP_WARNING = re.compile(r"^\*\*\*\s+(?:\S+:\s+)?(?:user )?[Ww]arning:")
# "***   Warning: increasing stack size", "***   user warning: ...", and the same raised inside a PARI
# function, which gp prefixes with the function name: "*** _^s: Warning: increasing stack size to ...".


class GPError(Exception):
    """gp failed, produced no sentinel, or reported an error."""


class GPTimeout(GPError):
    """gp exceeded the stage's wall-clock cap, or the run's deadline arrived first."""


_DEADLINE: contextvars.ContextVar = contextvars.ContextVar("ec_gp_deadline", default=None)


@contextlib.contextmanager
def deadline(at: float | None):
    """Clip every run_gp timeout inside the block to the time left until `at` (a time.monotonic()
    instant); None leaves the per-stage caps alone. The board budget is a hard limit this way: one
    slow curve cannot carry a run past it by exhausting each stage's own cap in turn."""
    token = _DEADLINE.set(at)
    try:
        yield
    finally:
        _DEADLINE.reset(token)


class GPVersionMismatch(GPError):
    pass


def repo_root() -> pathlib.Path:
    return pathlib.Path(__file__).resolve().parent.parent


def find_gp() -> str:
    """EC_GP_BIN, then the pinned build under tools/pari, then PATH."""
    env = os.environ.get("EC_GP_BIN")
    if env:
        return env
    built = repo_root() / "tools" / "pari" / "bin" / "gp"
    if built.is_file() and os.access(built, os.X_OK):
        return str(built)
    found = shutil.which("gp")
    if found:
        return found
    raise GPError("GP_UNAVAILABLE: no gp binary (set EC_GP_BIN or run ./setup.sh)")


def _env_for(gp_bin: str, scratch: str | None = None) -> dict:
    """Minimal, deterministic environment for gp. A source-built gp (tools/pari) may be
    dynamically linked against the libpari beside it; point the loader there."""
    env = {"PATH": os.environ.get("PATH", ""), "LC_ALL": "C"}
    if scratch:
        env.update({"TMPDIR": scratch, "HOME": scratch})
    if os.environ.get("GP_DATA_DIR"):
        env["GP_DATA_DIR"] = os.environ["GP_DATA_DIR"]
    libdir = pathlib.Path(gp_bin).resolve().parent.parent / "lib"
    if libdir.is_dir():
        env["LD_LIBRARY_PATH"] = str(libdir)
        env["DYLD_LIBRARY_PATH"] = str(libdir)
    return env


def gp_version(gp_bin: str) -> str:
    try:
        out = subprocess.run([gp_bin, "--version-short"], capture_output=True, text=True, timeout=30,
                             env=_env_for(gp_bin))
    except (OSError, subprocess.TimeoutExpired) as e:
        raise GPError(f"GP_UNAVAILABLE: cannot run {gp_bin}: {e}")
    return (out.stdout or out.stderr).strip()


def require_pinned(gp_bin: str) -> str:
    v = gp_version(gp_bin)
    if v != PINNED_VERSION:
        raise GPVersionMismatch(f"GP_VERSION_MISMATCH: gp reports {v!r}, benchmark pins {PINNED_VERSION}")
    return v


def gp_int(n: int) -> str:
    """Render a Python int for GP; refuses anything that is not an int."""
    if not isinstance(n, int) or isinstance(n, bool):
        raise TypeError("gp_int expects an int")
    s = str(n)
    if not _INT_RE.match(s):
        raise TypeError("gp_int: not an integer literal")
    return s


def gp_rat(num: int, den: int) -> str:
    return f"{gp_int(num)}/{gp_int(den)}" if den != 1 else gp_int(num)


def gp_vec(items) -> str:
    return "[" + ", ".join(items) + "]"


def gp_point(x, y) -> str:
    """A rational point [x, y] (Fractions) for GP."""
    return gp_vec([gp_rat(x.numerator, x.denominator), gp_rat(y.numerator, y.denominator)])


def gp_command(gp_bin: str, script: pathlib.Path) -> list:
    """The gp invocation for one stage. Stack sizes go on the command line (a source-built gp
    abandons the rest of a script file after an in-script parisize change). recover=0 makes gp
    abort the script on any error with a non-zero exit and no sentinel, so an error fails the
    stage structurally; the stderr classifier in run_gp is the second line and names the error."""
    return [gp_bin, "-q", "-f", "--default", f"parisize={PARISIZE}", "--default", f"parisizemax={PARISIZEMAX}",
            "--default", "recover=0", str(script)]


def run_gp(body: str, timeout: float, precision: int = DEFAULT_PRECISION, gp_bin: str | None = None) -> dict:
    """Run `body` (GP statements, each on its own line) and return the KEY=VALUE map."""
    gp_bin = gp_bin or find_gp()
    at = _DEADLINE.get()
    clipped = False
    if at is not None:
        remaining = at - time.monotonic()
        if remaining <= 0:
            raise GPTimeout("board budget exhausted before this gp stage started")
        clipped = remaining < timeout
        timeout = min(timeout, remaining)
    script = (
        f"default(realprecision, {int(precision)});\n"
        + body.rstrip("\n") + "\n"
        + f'print("{SENTINEL}=1");\n'
        "quit;\n"
    )
    scratch = tempfile.mkdtemp(prefix="ecgp-", dir=os.environ.get("TMPDIR") or None)
    path = pathlib.Path(scratch) / "stage.gp"
    try:
        path.write_text(script, encoding="ascii")
        try:
            env = _env_for(gp_bin, scratch)
            proc = subprocess.run(gp_command(gp_bin, path), capture_output=True, text=True, timeout=timeout, env=env, cwd=scratch)
        except subprocess.TimeoutExpired:
            if clipped:                        # the board ran out of time, not this curve's stage
                raise GPTimeout(f"board budget exhausted during this gp stage ({timeout:.0f}s were left)")
            raise GPTimeout(f"gp exceeded {timeout:.0f}s")
        except OSError as e:
            raise GPError(f"GP_UNAVAILABLE: {e}")
    finally:
        shutil.rmtree(scratch, ignore_errors=True)
    # With recover=0 an error aborts gp (exit 1, no sentinel), so the sentinel check below fails
    # closed on its own; this pass names the error and also covers a gp that kept going: every `***`
    # diagnostic on stderr that is not a warning (stack growth, user warnings) fails the stage.
    errors = [l.strip() for l in proc.stderr.splitlines()
              if l.strip().startswith("***") and not _GP_WARNING.match(l.strip())]
    if errors:
        raise GPError(f"GP_ERROR: gp reported an error: {' | '.join(errors)[:600]}")
    out: dict = {}
    for line in proc.stdout.splitlines():
        if "=" in line:
            k, v = line.split("=", 1)
            out[k.strip()] = v.strip()
    if out.get(SENTINEL) != "1":
        err = "\n".join(l for l in (proc.stderr + "\n" + proc.stdout).splitlines() if "***" in l)
        tail = (proc.stdout[-300:] + " | stderr: " + proc.stderr[-300:]).replace("\n", " ")
        raise GPError(f"GP_ERROR: gp did not complete (exit {proc.returncode}): {err[:600] or tail}")
    return out
