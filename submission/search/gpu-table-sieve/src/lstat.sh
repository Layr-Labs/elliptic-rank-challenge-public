#!/bin/bash
# lstat.sh - one-screen status of the lottery restart (dirs runs/camp_*_r*): sieves alive, GPU, candidates swept since the restart, ranked rank histogram, new events.
# usage (Spark, from ~/ecwork_src/kub):  ./lstat.sh [SUFFIX_REGEX=_r[0-9]+$]
cd "$(dirname "$0")"
echo "$(date -u +%H:%MZ) sieves=$(ps -eo args | grep -c '[k]ubg_q[0-9]* sieve') gpu=$(nvidia-smi --query-gpu=utilization.gpu --format=csv,noheader | tr -d ' ') load=$(cut -d' ' -f1 /proc/loadavg) rankers=$(ps -eo args | grep -c '[k]ub_rank.gp') inc_loops=$(ps -eo args | grep -c '[i]nc_loop.sh')"
python3 - "${1:-_r[0-9]+$|camp_(44|66|67|5)x}" <<'PY'
import glob, re, collections, sys
pat = re.compile(sys.argv[1])
tot = 0.0; done = 0; n = 0; alive_lines = 0
cnt = collections.Counter(); r9 = {}; hi = []
rows = {}
for d in sorted(glob.glob("runs/camp_*")):
    if not pat.search(d): continue
    n += 1
    try: L = open(d + "/sieve.log", errors="ignore").read().split("\n")
    except OSError: continue
    last = [l for l in L if l.startswith(("q<", "p<"))]
    dn = [l for l in L if l.startswith("done:")]
    if dn:
        done += 1
        m = re.search(r"candidates (\d+)", dn[-1])
        if m: tot += float(m.group(1))
    elif last:
        m = re.search(r"cand ([0-9.e+]+)", last[-1])
        if m: tot += float(m.group(1))
    for f in glob.glob(d + "/inc/ranked_*.txt") + glob.glob(d + "/ranked.txt"):
        for l in open(f, errors="ignore"):
            t = l.split()
            if len(t) < 10 or not re.fullmatch(r"-?\d+", t[0]) or not re.fullmatch(r"-?\d+", t[6]): continue
            try: rows[(t[0], t[1])] = (int(t[6]), int(t[7]), float(t[5]), d)
            except ValueError: pass
err = [k for k, (lb, ub, H, d) in rows.items() if lb == -1]
for k, (lb, ub, H, d) in rows.items():
    if lb >= 6: cnt[lb] += 1
    if lb == 9: r9[k] = round(H, 3)
    if lb >= 10 or ub >= 10: hi.append((k, lb, ub, H, d))
print("dirs=%d done=%d swept=%.3g cand | ranked>=6: %s | r9 H: %s | r10+/ub10: %s | rank errors/timeouts (lb=-1): %d" % (n, done, tot, dict(sorted(cnt.items())), sorted(r9.values()), hi, len(err)))
# per progression class (camp_44x{a|b}{T}n{N}s{j}_{HLO}_{HHI}): candidates swept, dirs done/total, curves with rank>=7 by rank
cls = collections.OrderedDict()
for d in sorted(glob.glob("runs/camp_*x*")):
    m = re.match(r"runs/camp_(\d+)x[ab]?(\d)n(\d)(?:s\d+)?_([0-9.]+)_([0-9.]+)", d)
    if not m or not pat.search(d): continue
    key = "F%s t%sn%s[%s,%s)" % m.groups()
    c = cls.setdefault(key, [0.0, 0, 0, collections.Counter()])
    try: L = open(d + "/sieve.log", errors="ignore").read().split("\n")
    except OSError: continue
    c[2] += 1
    dn = [l for l in L if l.startswith("done:")]
    last = [l for l in L if l.startswith(("q<", "p<"))]
    if dn:
        c[1] += 1; mm = re.search(r"candidates (\d+)", dn[-1]); c[0] += float(mm.group(1)) if mm else 0
    elif last:
        mm = re.search(r"cand ([0-9.e+]+)", last[-1]); c[0] += float(mm.group(1)) if mm else 0
    for (p_, q_), (lb, ub, H, dd) in rows.items():
        if dd == d and lb >= 7: c[3][lb] += 1
for key, (cd, dn_, nd, rc) in cls.items():
    print("  %-24s cand %.2e dirs %d/%d done | r7 %d r8 %d r9 %d r10+ %d" % (key, cd, dn_, nd, rc.get(7, 0), rc.get(8, 0), rc.get(9, 0), sum(v for k, v in rc.items() if k >= 10)))
PY
python3 ev_new.py 2>&1 | tail -n 6
