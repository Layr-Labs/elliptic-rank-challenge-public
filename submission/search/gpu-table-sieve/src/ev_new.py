#!/usr/bin/env python3
"""ev_new.py - print only NEW high-rank events of the ranked files under ~/ecwork_src/kub/runs (state in events_seen.txt).
Events: any curve with lower or upper bound >= 10; or lower bound 9 below the current cell holder (Legendre 4x2/9 19.0165, Z/6 18.446, Z/5 20.47, Z/4 12.2).
usage: ev_new.py [--init]   (--init: mark everything present as seen without printing)"""
import glob, os, re, sys
ROOT = os.path.expanduser("~/ecwork_src/kub/runs")
SEEN = os.path.expanduser("~/ecwork_src/kub/events_seen.txt")
seen = set(open(SEEN).read().split("\n")) if os.path.exists(SEEN) else set()
new = []
for d in sorted(glob.glob(ROOT + "/camp_*") + glob.glob(ROOT + "/z4_*")):
    z4 = os.path.basename(d).startswith("z4_")
    fam = "Z4" if z4 else os.path.basename(d).split("_")[1]
    for f in glob.glob(d + "/inc/ranked_*.txt") + glob.glob(d + "/ranked.txt"):
        for l in open(f, errors="ignore"):
            t = l.split()
            if len(t) < 10 or not re.fullmatch(r"-?\d+", t[0]) or not re.fullmatch(r"-?\d+", t[6]) or not re.fullmatch(r"-?\d+", t[7]):
                continue
            lb, ub = int(t[6]), int(t[7])
            try:
                H = float(t[5])
            except ValueError:
                continue
            # thresholds: an event = rank >= 10 (breakthrough cell), or rank 9 below the current cell holder (Legendre 19.02, Z/6 18.445, Z/5 20.47, Z/4 12.2)
            if z4: h9 = 12.2
            elif fam in ("66", "67") or fam.startswith("66") or fam.startswith("67"): h9 = 18.446
            elif fam.startswith("5"): h9 = 20.47
            else: h9 = 19.0165
            hit = lb >= 10 or ub >= 10 or (lb == 9 and H < h9)
            if not hit:
                continue
            key = "%s %s %s" % ("C" if not z4 else "Z4", t[0], t[1])     # same curve found by several campaigns = one event
            if key in seen:
                continue
            seen.add(key)
            new.append("%s H=%.4f lb=%d ub=%d  (%s)" % (key, H, lb, ub, os.path.basename(d)))
    for e in glob.glob(d + "/found/*.ec"):
        if re.match(r"^[0-9]+-9-", os.path.basename(e)):
            continue                      # rank-9 files are routine (the ranked-line events apply the cell-holder heights); rank >= 10 files are reported
        key = "EC " + os.path.basename(e)          # file name = family-npts-p_q: the same curve in several campaigns is one event
        if key not in seen:
            seen.add(key); new.append(key + "  (%s)" % os.path.basename(d))
with open(SEEN, "w") as fh:
    fh.write("\n".join(sorted(seen)) + "\n")
if "--init" not in sys.argv:
    for n in new:
        print(n)
    if new:
        print("NEW EVENTS", len(new))
