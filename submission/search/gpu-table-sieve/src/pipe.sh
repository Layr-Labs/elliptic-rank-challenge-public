#!/bin/bash
# pipe.sh - sieve (GPU) -> rank stage (PARI on the CPU cores) -> summary -> verifier screen.   Runs on the Spark, from ~/ecwork_src/kub.
# usage: ./pipe.sh FAM HMAX QLO QHI THR1 THR2 THR3 OUTDIR [MINRANK=99] [NGP=16] [SIEVE=kubg]
#   FAM 5 | 6 | 42 | 8 | 10 | 12 | 26 | 28 (see fams.txt);  thresholds are on S(300), S(2000), S(10^4);  every output file lives in OUTDIR
#   MINRANK: curves whose ellrank lower bound reaches it are saturated and written to OUTDIR/found/*.ec, then screened with ec_verify
#   env KUB_HMIN=<h>: shell lower bound;  env PIPE_ASYNC_RANK=1: run the rank stage in the background (so the next sieve can start at once)
# bounded by construction: one sieve run, one rank pass; no loops, no cron.
set -u
FAM=$1; HMAX=$2; QLO=$3; QHI=$4; T1=$5; T2=$6; T3=$7; OUT=$8
MINRANK=${9:-99}; NGP=${10:-16}; SIEVE=${11:-kubg}
HERE=$(cd "$(dirname "$0")" && pwd)
GP=$HOME/ecrank/tools/pari/bin/gp
export KUB_TABDIR=$HERE/tabs
export KUB_FAMS=$HERE/fams.txt
[ -f "$HERE/pipe.conf" ] && . "$HERE/pipe.conf"        # runtime-tunable settings (PIPE_MINS: rank only candidates with S(10^4) >= PIPE_MINS, park the rest in OUT/lo_backlog.txt)
mkdir -p "$OUT/found"
OUT=$(cd "$OUT" && pwd)          # absolute path: the verifier screen below runs from ~/ecrank
rm -f "$OUT/rank_done" "$OUT/FOUND.txt"
t0=$(date +%s)
if [ "${PIPE_SKIP_SIEVE:-0}" = "1" ]; then echo "sieve skipped: using the existing $OUT/cand.txt"      # env PIPE_SKIP_SIEVE=1: rank-only run on an existing candidate list
elif [[ "$SIEVE" == kubg* ]]; then "$HERE/$SIEVE" sieve "$FAM" "$HMAX" "$QLO" "$QHI" "$T1" "$T2" "$T3" "${PIPE_HOSTTHR:-6}" > "$OUT/cand.txt" 2> "$OUT/sieve.log"      # any GPU sieve binary named kubg* (e.g. kubg_g = generalized a4/a6 families)
else "$HERE/kub" sieve "$FAM" "$HMAX" "$QLO" "$QHI" "$T1" "$T2" "$T3" 20 > "$OUT/cand.txt" 2> "$OUT/sieve.log"; fi
t1=$(date +%s)
n=$(wc -l < "$OUT/cand.txt")
echo "sieve: $(tail -1 "$OUT/sieve.log" 2>/dev/null) ($((t1-t0)) s wall, $n candidates to rank)"

rank_stage() {
  rm -f "$OUT"/part_*
  for _w in $(seq 1 1440); do                                                    # an inc_rank.sh call (incremental ranking while the sieve ran) may still be running: wait for it (up to 2 h)
    [ -f "$OUT/inc_lock" ] || break
    _lp=$(cat "$OUT/inc_lock" 2>/dev/null); if [ -n "$_lp" ] && ! kill -0 "$_lp" 2>/dev/null; then break; fi
    sleep 5
  done
  inc=0; [ -f "$OUT/inc_done" ] && inc=$(cat "$OUT/inc_done")                    # lines of cand.txt already ranked by inc_rank.sh
  tail -n +$((inc+1)) "$OUT/cand.txt" > "$OUT/cand.rest"
  if [ "${PIPE_MINS:-0}" != "0" ]; then
    awk -v h="$PIPE_MINS" '$5 <  h' "$OUT/cand.rest" >> "$OUT/lo_backlog.txt"
    awk -v h="$PIPE_MINS" '$5 >= h' "$OUT/cand.rest" > "$OUT/cand.rest2"; mv "$OUT/cand.rest2" "$OUT/cand.rest"
  fi
  n=$(wc -l < "$OUT/cand.rest")
  if [ "$n" -gt 0 ]; then
    if [ "$(head -1 "$OUT/cand.rest" | wc -w)" -ge 7 ]; then sort -k7 -n -r "$OUT/cand.rest" > "$OUT/cand.sorted"; else sort -k5 -n -r "$OUT/cand.rest" > "$OUT/cand.sorted"; fi          # best first: S(1e5) if present, else S(10^4)
    split -n r/"$NGP" -d "$OUT/cand.sorted" "$OUT/part_"           # round robin: every process sees the best candidates first
    for p in "$OUT"/part_??; do
      ( KUB_FAM=$FAM KUB_IN=$p KUB_MINRANK=$MINRANK EC_OUT=$OUT/found KUB_TMAX=${KUB_TMAX:-600} nice -n 5 "$GP" -q "$HERE/kub_rank.gp" > "$p.out" 2> "$p.err" ) &
    done
    wait
  fi
  t2=$(date +%s)
  { cat "$OUT"/inc/ranked_*.txt 2>/dev/null; cat "$OUT"/part_??.out 2>/dev/null | grep -v '^ \*\*\*'; } > "$OUT/ranked.txt"
  echo "rank stage: $((t2-t1)) s wall"
  python3 - "$OUT/ranked.txt" <<'PYEOF'
import sys, collections
rows = []
for line in open(sys.argv[1]):
    t = line.split()
    if len(t) < 10 or t[2] == "torsion-mismatch": continue
    try: rows.append((int(t[0]), int(t[1]), float(t[2]), float(t[3]), float(t[4]), float(t[5]), int(t[6]), int(t[7])))
    except ValueError: pass
c = collections.Counter(r[6] for r in rows)
print("ranked", len(rows), "curves; lower-bound rank histogram:", sorted(c.items()))
for k in sorted(c):
    if k >= 5:
        hs = sorted(r[5] for r in rows if r[6] == k)
        print("  rank >= %d: n=%d  H min %.3f  median %.3f" % (k, c[k], hs[0], hs[len(hs)//2]))
PYEOF
  if ls "$OUT"/found/*.ec >/dev/null 2>&1; then
    ( cd "$HOME/ecrank" && export PATH="$HOME/ecrank/tools/pari/bin:$PATH"
      for f in "$OUT"/found/*.ec; do
        s0=$(date +%s.%N); r=$(.venv-bench/bin/python -m ec_verify "$f" 2>&1 | tail -2 | tr '\n' ' '); s1=$(date +%s.%N)
        printf "%s  [%.1f s]  %s\n" "$(basename "$f")" "$(echo "$s1 - $s0" | bc)" "$r"
      done ) > "$OUT/verify.log" 2>&1
    echo "verifier screen: $(wc -l < "$OUT/verify.log") curves, see $OUT/verify.log"
    cp "$OUT/verify.log" "$OUT/FOUND.txt"
  fi
  touch "$OUT/rank_done"
}

if [ "${PIPE_ASYNC_RANK:-0}" = "1" ]; then
  ( rank_stage ) > "$OUT/rank_report.txt" 2>&1 < /dev/null &     # the subshell itself must not hold the caller's pipe open (campaign.sh pipes us into tail)
  echo "rank stage started in the background (report: $OUT/rank_report.txt)"
else
  rank_stage
fi
