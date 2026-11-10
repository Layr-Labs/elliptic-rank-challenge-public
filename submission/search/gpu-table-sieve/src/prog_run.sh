#!/bin/bash
# prog_run.sh - launch one PROGRESSION-CLASS cell of the Legendre (FAM 44) rank-10 hunt: outer value with EXACTLY T of the primes 7..47, inner value divisible by at least N
# further primes of that set (KUB_PROG=N, N<=3), naive height shell [HLO,HHI), both roles A and B, S slices (KUB_QMOD) per role.   Runs on the Spark, from ~/ecwork_src/kub.
# usage: [PROG_FAM=44|5|66|67] ./prog_run.sh T N HLO HHI [SLICES=2] [QHI] [KUBG=kubg_q9] [GAP=2]      (FAM 44 / 5: roles A and B; FAM 66 / 67 (Z/6, outer p resp. p+q): one pass, no role)
# Why: measured on our own rank 6-9 finds, a requirement "inner has >= 2 further small bad primes" keeps ~20-45 % of the high-rank curves in ~3-5 % of the candidates (x5-9 per candidate; test cell
# t=5 n=2 [21.5,22.5): 3.3e13 candidates in 17 min -> 10 rank-8, plain sweep needs 8e14 candidates for ~46), and with items of 1e7+ candidates (outer q >= 1e9) the per-item overhead is nil.
# Output dirs runs/camp_44x<a|b><T>n<N>s<j>_<HLO>_<HHI> (tag x = progression class), logs runs/campaign_<tag>.log.  Finite: every process is one pipe.sh run, no loops, no cron.
set -u
cd "$(dirname "$0")"
T=$1; N=$2; HLO=$3; HHI=$4; S=${5:-2}; KUBG=${7:-kubg_q9}; GAP=${8:-2}
QHI=${6:-$(python3 -c "import math; print(int(2*math.exp($HHI-0.69)))")}
FAM=${PROG_FAM:-44}
if [ "$FAM" = 44 ] || [ "$FAM" = 5 ]; then ROLES="A B"; else ROLES="-"; fi
for role in $ROLES; do
  r=$(echo $role | tr A-Z a-z); roleenv="KUB_ROLE=$role"; if [ "$role" = "-" ]; then r=""; roleenv=""; fi
  for ((j = 0; j < S; j++)); do
    tag=x${r}${T}n${N}s${j}; [ "$S" = 1 ] && tag=x${r}${T}n${N}
    qm=""; [ "$S" -gt 1 ] && qm="KUB_QMOD=$S:$j"
    env CAMP_KUBG=$KUBG CAMP_TAG=$tag CAMP_QHI=$QHI KUB_QRICH=$T KUB_QRICHMAX=$T KUB_PROG=$N $roleenv $qm KUB_THR4=55 KUB_TMAX=900 \
      setsid nohup ./campaign.sh $FAM 9.5 24 34 9 $HLO $HHI > runs/campaign_${FAM}${tag}_${HLO}_${HHI}.log 2>&1 < /dev/null &
    echo "launched $tag [$HLO,$HHI) QHI=$QHI"
    sleep "$GAP"
  done
done
