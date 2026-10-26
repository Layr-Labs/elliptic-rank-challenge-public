#!/bin/bash
# usage: z4stream.sh outdir S10kmin ; ranks finished chunks; exits when all jobs ranked
OUT=$1; TH=$2
while true; do
  NJ=$(wc -l < $OUT/jobs.txt); ND=0
  for d in $OUT/c.*.txt.done $OUT/c.*[0-9].done; do
    [ -e "$d" ] || continue; case $d in *.r.done) continue;; esac
    b=${d%.done}; ND=$((ND+1))
    [ -f $b.r.done ] && continue
    awk -v t=$TH '$7>=t' $b.txt > $b.sel
    if [ -s $b.sel ]; then printf 'read("rank4.gp")\nrunrank4("%s","%s",15.39)\n' $b.sel $b.r | gp -q --default parisizemax=1000000000 >/dev/null 2>&1; fi
    touch $b.r.done
  done
  [ $ND -ge $NJ ] && break
  sleep 120
done
echo z4stream-done $OUT
