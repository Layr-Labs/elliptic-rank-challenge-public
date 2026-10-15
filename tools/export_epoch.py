#!/usr/bin/env python3
"""Write a submission in Epoch AI's FrontierMath 'Elliptic Curves over Q of Large Rank' JSON format:
{"curve": [a1..a6 as strings], "x_coords": [x as integer or fraction strings]}.

    python tools/export_epoch.py submission/board/<file>.ec > <file>.epoch.json
"""

import json
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from ec_verify.parse import parse  # noqa: E402


def main() -> int:
    s = parse(pathlib.Path(sys.argv[1]).read_text(encoding="ascii"))
    xs = [str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}" for x, _ in s.points]
    json.dump({"curve": [str(a) for a in s.ainvs], "x_coords": xs}, sys.stdout)
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
