#!/usr/bin/env python3
"""Owner tool: turn a ranked line of famgate.gp / z8_gate.gp into a .ec file (repo/README.md format).
usage: mkec.py <ranked-file> <fam> <p> <q> <out.ec> [comment...]"""
import sys, re
rf, fam, p, q, out = sys.argv[1:6]
comment = " ".join(sys.argv[6:])
for line in open(rf):
    t = line.split()
    if len(t) > 3 and t[0] == fam and t[1] == p and t[2] == q:
        m = re.search(r" ainv (\[[^\]]*\]) pts (\[.*\])\s*$", line)
        ainv = [x.strip() for x in m.group(1).strip("[]").split(",")]
        pts = re.findall(r"\[\s*([-0-9/]+)\s*,\s*([-0-9/]+)\s*\]", m.group(2))
        with open(out, "w") as f:
            if comment: f.write("# %s\n" % comment)
            f.write("E " + " ".join(ainv) + "\n")
            for x, y in pts: f.write("P %s %s\n" % (x, y))
            f.write("#CLAIM rank=%d\n" % len(pts))
        print("wrote", out, "with", len(pts), "points")
        break
else:
    sys.exit("line not found")
