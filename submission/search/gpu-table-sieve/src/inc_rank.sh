#!/bin/bash
# inc_rank.sh - rank the candidates that a RUNNING sieve has appended to DIR/cand.txt since the last call (so a hit shows up within minutes, not at the end of the shell).
# usage: ./inc_rank.sh DIR FAM MINRANK [NGP=18]      (DIR = runs/camp_<FAM>_<Hlo>_<Hhi>[_x]; run from ~/ecwork_src/kub; env KUB_DSMIN default 8, KUB_WSIGN, KUB_TMAX, INC_HI)
# Progress is kept in DIR/inc_done (number of cand.txt lines already ranked); chunk results in DIR/inc/; pipe.sh's final rank stage ranks only the lines after inc_done.
# INC_HI=<s> (default 0 = off): within a chunk the lines with S(10^4) >= s are ranked first (phase A, FOUND.txt written right after it), the rest afterwards (phase B),
# or, with INC_SKIPLO=1, parked in DIR/lo_backlog.txt.
# A curve with lower-bound rank >= MINRANK lands in DIR/found/*.ec, is screened with ec_verify (novelty ON) and DIR/FOUND.txt is (re)written.
set -u
DIR=$(cd "$1" && pwd); FAM=$2; MINRANK=$3; NGP=${4:-18}
HERE=$(cd "$(dirname "$0")" && pwd)
GP=$HOME/ecrank/tools/pari/bin/gp
export KUB_FAMS=$HERE/fams.txt KUB_DSMIN=${KUB_DSMIN:-8}
HI=${INC_HI:-0}
mkdir -p "$DIR/inc" "$DIR/found"
if [ -f "$DIR/inc_lock" ]; then _lp=$(cat "$DIR/inc_lock" 2>/dev/null); if [ -z "$_lp" ] || kill -0 "$_lp" 2>/dev/null; then echo "inc_rank: another call is running (inc_lock)"; exit 0; fi; fi
echo $$ > "$DIR/inc_lock"; trap 'rm -f "$DIR/inc_lock"' EXIT
done_n=$(cat "$DIR/inc_done" 2>/dev/null || echo 0)
total=$(wc -l < "$DIR/cand.txt")                 # complete lines only
if [ "$total" -le "$done_n" ]; then echo "inc_rank: nothing new ($total lines)"; exit 0; fi
tag=$(printf "%09d" "$done_n")
tail -n +$((done_n+1)) "$DIR/cand.txt" | head -n $((total-done_n)) > "$DIR/inc/chunk_$tag.txt"

verify_found() {
  if ls "$DIR"/found/*.ec >/dev/null 2>&1; then
    ( cd "$HOME/ecrank" && export PATH="$HOME/ecrank/tools/pari/bin:$PATH"
      for f in "$DIR"/found/*.ec; do
        r=$(.venv-bench/bin/python -m ec_verify "$f" 2>&1 | tail -2 | tr '\n' ' ')
        printf "%s  %s\n" "$(basename "$f")" "$r"
      done ) > "$DIR/verify.log" 2>&1
    cp "$DIR/verify.log" "$DIR/FOUND.txt"
    echo "inc_rank: !! FOUND curves with rank >= $MINRANK, see $DIR/FOUND.txt"; cut -c1-300 "$DIR/FOUND.txt"
  fi
}

rank_phase() {   # rank_phase INFILE NAME : sorted by S(10^4), round robin over NGP gp processes
  local in=$1 name=$2
  [ -s "$in" ] || return 0
  if [ "$(head -1 "$in" | wc -w)" -ge 7 ]; then sort -k7 -n -r "$in" > "$in.sorted"; else sort -k5 -n -r "$in" > "$in.sorted"; fi     # best candidates first: S(1e5) if the sieve printed it, else S(1e4)
  split -n r/"$NGP" -d "$in.sorted" "$DIR/inc/part_${name}_"
  for p in "$DIR"/inc/part_${name}_??; do
    [ -s "$p" ] || continue
    ( KUB_FAM=$FAM KUB_IN=$p KUB_MINRANK=$MINRANK EC_OUT=$DIR/found KUB_TMAX=${KUB_TMAX:-600} nice -n 5 "$GP" -q "$HERE/kub_rank.gp" > "$p.out" 2> "$p.err" ) &
  done
  wait
  cat "$DIR"/inc/part_${name}_??.out 2>/dev/null | grep -v '^ \*\*\*' >> "$DIR/inc/ranked_$tag.txt"
  rm -f "$DIR"/inc/part_${name}_?? "$DIR"/inc/part_${name}_??.err "$in.sorted"
}

t0=$(date +%s)
: > "$DIR/inc/ranked_$tag.txt"
if awk -v h="$HI" 'BEGIN{exit !(h>0)}'; then
  awk -v h="$HI" '$5 >= h' "$DIR/inc/chunk_$tag.txt" > "$DIR/inc/chunk_${tag}_hi.txt"
  awk -v h="$HI" '$5 < h'  "$DIR/inc/chunk_$tag.txt" > "$DIR/inc/chunk_${tag}_lo.txt"
  rank_phase "$DIR/inc/chunk_${tag}_hi.txt" "${tag}h"
  verify_found
  if [ "${INC_SKIPLO:-0}" = "1" ]; then cat "$DIR/inc/chunk_${tag}_lo.txt" >> "$DIR/lo_backlog.txt"          # low-S lines are parked, ranked later if there is time (sorted by S(10^4))
  else rank_phase "$DIR/inc/chunk_${tag}_lo.txt" "${tag}l"; fi
else
  rank_phase "$DIR/inc/chunk_$tag.txt" "$tag"
fi
echo "$total" > "$DIR/inc_done"
awk -v a="$((done_n+1))" -v b="$total" -v w="$(( $(date +%s) - t0 ))" '{c[$7]++} END{printf "inc_rank: lines %s-%s ranked in %s s; lower-bound histogram (rank:count, -2 = screened out by dS, -3 = parity, -1 = ellrank error or timeout): ", a, b, w; for (k in c) printf "%s:%d ", k, c[k]; printf "\n"}' "$DIR/inc/ranked_$tag.txt"
verify_found
# curves whose ellrank failed (error/timeout): list them for a second attempt
awk '$7 == -1' "$DIR/inc/ranked_$tag.txt" >> "$DIR/inc/errors.txt"
