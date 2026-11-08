#!/bin/bash
# ladder.sh FAM [T1 T2 T3] - bounded density ladder: five samples of the family at increasing integral-model height, identical cuts and the dS>=8 screen,
# everything that passes is ranked (MINRANK 9 -> only a rank-9 curve is written as .ec).  q-subsampling KUB_QMOD=m:1 thins the large bins.
# Env: LADDER_QM="m1 m2 m3 m4 m5" (default "1 7 47 163 641", sized for ~2.6e12 candidates per bin of the Z/5 family), LADDER_KUBG=<sieve binary> (default kubg),
#      LADDER_TAG=<suffix for the run dir>, LADDER_H="h0 h1 h2 h3 h4 h5" (bin edges, default "14.4 15.0 15.8 16.8 17.5 18.2").
# Output: runs/ladder_<FAM><TAG>/L1..L5 and ladder.log lines (rank histograms per bin).  Own code; used to measure rank-7/8 densities per family before a long campaign.
cd ~/ecwork_src/kub
FAM=${1:?usage: ladder.sh FAM [T1 T2 T3]}; T1=${2:-9.5}; T2=${3:-24}; T3=${4:-38}
export KUB_DSMIN=${KUB_DSMIN:-8}
read -r -a QMS <<< "${LADDER_QM:-1 7 47 163 641}"
read -r -a HE <<< "${LADDER_H:-14.4 15.0 15.8 16.8 17.5 18.2}"
KBIN=${LADDER_KUBG:-kubg}
OUTD=runs/ladder_${FAM}${LADDER_TAG:-}; mkdir -p $OUTD
run() { # name hmin hmax qmod
  local n=$1 lo=$2 hi=$3 qm=$4
  echo "== FAM $FAM bin $n H in [$lo,$hi) QMOD=$qm start $(date -u +%T)"
  if [ "$qm" = "1" ]; then KUB_HMIN=$lo ./pipe.sh $FAM $hi 1 1000000000 $T1 $T2 $T3 $OUTD/$n 9 18 $KBIN 2>&1 | tail -6
  else KUB_HMIN=$lo KUB_QMOD=$qm:1 ./pipe.sh $FAM $hi 1 1000000000 $T1 $T2 $T3 $OUTD/$n 9 18 $KBIN 2>&1 | tail -6; fi
}
for i in 0 1 2 3 4; do run L$((i+1)) ${HE[$i]} ${HE[$((i+1))]} ${QMS[$i]}; done
echo "== ladder FAM $FAM done $(date -u +%T)"
