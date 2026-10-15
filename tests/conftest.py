import pathlib
import shutil
import sys

import pytest

ROOT = pathlib.Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from ec_verify import gp as gpmod  # noqa: E402


def _gp_ok() -> bool:
    try:
        gpmod.require_pinned(gpmod.find_gp())
        return True
    except gpmod.GPError:
        return False


needs_gp = pytest.mark.skipif(not _gp_ok(), reason=f"pinned gp {gpmod.PINNED_VERSION} not available (run ./setup.sh)")
