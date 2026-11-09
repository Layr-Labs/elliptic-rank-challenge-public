#!/bin/bash
# campaign.sh - bounded rank-9 hunt: increasing-H shells for one family, each shell = sieve (GPU) -> rank -> ec_verify screen.
# usage: ./campaign.sh FAM THR1 THR2 THR3 MINRANK H0 H1 H2 ...      (shell k covers H in [H(k-1), H(k)); the first list entry is the lower bound of shell 1)
#   e.g.  ./campaign.sh 6 10 26 42 9 15.3 16.5 17.3 17.9 18.45        (Z/6, shells [15.3,16.5) [16.5,17.3) [17.3,17.9) [17.9,18.45))
# The rank stage of shell k runs in the background while the GPU already sieves shell k+1.
# Finite by construction (no loops beyond the argument list), no cron, no daemons. Results: runs/camp_<FAM>_<Hlo>_<Hhi>/{cand.txt,ranked.txt,found/*.ec,verify.log}
set -u
FAM=$1; T1=$2; T2=$3; T3=$4; MINRANK=$5; shift 5
HERE=$(cd "$(dirname "$0")" && pwd)
cd "$HERE"
prev=$1; shift
for hi in "$@"; do
  out=runs/camp_${FAM}${CAMP_TAG:-}_${prev}_${hi}          # CAMP_TAG distinguishes parallel campaigns of one family (e.g. qa3 / qb3 for the outer-value filter roles)
  echo "== shell FAM=$FAM H in [$prev,$hi)  start $(date -u +%T)"
  KUB_HMIN=$prev PIPE_ASYNC_RANK=1 ./pipe.sh "$FAM" "$hi" 1 "${CAMP_QHI:-1000000000}" "$T1" "$T2" "$T3" "$out" "$MINRANK" 18 "${CAMP_KUBG:-kubg}" 2>&1 | tail -12
  if ls "$out"/found/*.ec >/dev/null 2>&1; then
    echo "!! FOUND curves with rank >= $MINRANK in shell [$prev,$hi):"; ls -la "$out"/found/; cat "$out"/verify.log 2>/dev/null
  fi
  prev=$hi
done
echo "== campaign done $(date -u +%T)"
