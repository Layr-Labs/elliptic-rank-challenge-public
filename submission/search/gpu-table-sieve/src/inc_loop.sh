#!/bin/bash
# inc_loop.sh - bounded driver for inc_rank.sh: while sieves of family FAM are running, rank what each has produced so far, every INTERVAL seconds, at most ITER times.
# usage: ./inc_loop.sh FAM MINRANK ITER INTERVAL      (run from ~/ecwork_src/kub; env KUB_DSMIN / KUB_WSIGN / KUB_TMAX are passed on to inc_rank.sh)
# Handles every running "kubg*/kub sieve FAM" process (several shells of one family may sieve at the same time): the output dir is taken from the pipe.sh that started it.
# Stops early when no sieve of FAM is running for two consecutive checks (the campaign's own rank stage ranks the rest) or when a FOUND.txt appears.
set -u
FAM=$1; MINRANK=$2; ITER=$3; INTERVAL=$4
cd "$(dirname "$0")"
idle=0
for i in $(seq 1 "$ITER"); do
  DIRS=$(ps -eo pid,ppid,args | grep "[k]ubg[_a-z0-9]* sieve $FAM " | awk '{print $2}' | while read KP; do ps -o args= -p "$KP" | awk '/pipe.sh/ {print $10}'; done | sort -u)
  if [ -z "$DIRS" ]; then idle=$((idle+1)); echo "$(date -u +%T) no running sieve for family $FAM (idle $idle)"; [ "$idle" -ge 2 ] && break; sleep "$INTERVAL"; continue; fi
  idle=0
  stop=0
  for OUT in $DIRS; do
    [ "${OUT#/}" = "$OUT" ] && OUT="$PWD/$OUT"
    echo "$(date -u +%T) iteration $i: $OUT"
    ./inc_rank.sh "$OUT" "$FAM" "$MINRANK" 18 2>&1 | tail -n 6
    if [ -f "$OUT/FOUND.txt" ]; then echo "$(date -u +%T) FOUND marker in $OUT"; stop=1; fi
  done
  [ "$stop" = "1" ] && break
  sleep "$INTERVAL"
done
echo "$(date -u +%T) inc_loop finished"
