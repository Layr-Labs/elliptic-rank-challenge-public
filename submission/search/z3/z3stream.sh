#!/bin/bash
# usage: z3stream.sh outdir S10kmin Hmax
OUT=$1; TH=$2; HM=$3
while true; do
  NJ=$(wc -l < $OUT/jobs.txt); ND=0
  for d in $OUT/c.*[0-9].done; do
    [ -e "$d" ] || continue; b=${d%.done}; ND=$((ND+1))
    [ -f $b.r.done ] && continue
    awk -v t=$TH '$6>=t' $b.txt > $b.sel
    if [ -s $b.sel ]; then printf 'read("rank3.gp")\nrunrank3("%s","%s",%s)\n' $b.sel $b.r $HM | gp -q --default parisizemax=1000000000 >/dev/null 2>&1; fi
    touch $b.r.done
  done
  [ $ND -ge $NJ ] && break
  sleep 120
done
echo z3stream-done $OUT
