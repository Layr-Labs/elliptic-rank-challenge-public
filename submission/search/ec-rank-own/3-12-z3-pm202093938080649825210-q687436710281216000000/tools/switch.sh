#!/bin/bash
# build-dside switch (20 Nov 2026): stock lotsieve7 -> lotsieve7v --vmul at the sweep's checkpoint, same outdir (survivors appended).
# usage: switch.sh STOCKPRE STOCKTAG STOCKCPUS     e.g. card1: enum/z3u12-o15r z3u12-o15r 0-23 ; card2: enum/z3u12-o15hi z3u12-o15hi 0-22
cd /root/ecr; O=out/z3u12-o15r; SP=$1; STAG=$2; SCPU=$3; V=enum/z3u12-o15v; VS=enum/z3u12-o15vs
SIEVE="--table search/lottery/tables/z3 --outdir $O --cut1 10 --cut2 28 --cut3 40 --hshell 14.5 17.1438 --batch 2e8 --stop /root/ecr/STOP --family z3w --no-watchdog-i-know --resume"
restart_stock() { rm -f STOP; setsid nohup taskset -c $SCPU ./bin/lotsieve7 $SIEVE --spans $SP --tag $STAG >> $O/run.log 2>&1 < /dev/null & echo "stock restarted"; }
[ -s $V.info ] || { echo "vmul enum $V not finished"; exit 1; }
pgrep -x lotsieve7v >/dev/null && { echo "lotsieve7v already running"; exit 1; }
ts=$(date -u +%FT%TZ); t0=$(date +%s)
touch STOP
for i in $(seq 1 120); do pgrep -x lotsieve7 >/dev/null || break; sleep 1; done
if pgrep -x lotsieve7 >/dev/null; then echo "stock still running after 120 s, abort"; rm -f STOP; exit 1; fi
cp $O/checkpoint $O/checkpoint.stock.$ts
if ! dside/vslice $V $VS $SP $O/checkpoint > dside/ck.new 2> dside/vslice.log; then echo "vslice failed:"; cat dside/vslice.log; restart_stock; exit 1; fi
cp dside/ck.new $O/checkpoint
rm -f STOP
echo "$ts build-dside: switched to dside/lotsieve7v --vmul (row rule v2(b)=1 -> 2|a, v3(b)=1 -> 3|a; validated on this card, dside/val) at stock checkpoint '$(cat $O/checkpoint.stock.$ts)'; spans $VS ($(sed -n 2p $VS.info)); $(cat dside/vslice.log)" >> $O/command.txt
setsid nohup taskset -c 0-22 dside/lotsieve7v $SIEVE --spans $VS --tag $STAG-vmul --vmul >> $O/run.log 2>&1 < /dev/null &
sleep 30; echo "downtime about $(( $(date +%s) - t0 - 30 )) s before relaunch"; pgrep -x lotsieve7v >/dev/null && echo "lotsieve7v running" || { echo "lotsieve7v NOT running"; tail -5 $O/run.log; }
tail -2 $O/progress.tsv | cut -f1-10; cat $O/checkpoint; cat dside/vslice.log
