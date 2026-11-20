#!/bin/bash
# build-dside (20 Nov 2026): when this card's 3/12 sweep (dside/lotsieve7v -> out/z3u12-o15r) exits, start the 3/13 omega=18 half
# (enum/z3u13-o18, --vmul) so the GPU does not idle (the lead's guard destroys a card after 30 min under 10 %).
# Decision on the status lines written AFTER this waiter started (run.log already holds old STOPPED lines from swaps):
#   new DONE -> chain;  new STOPPED or STOP file -> no chain (deliberate stop);  no new status line (crash) -> log CRASH and chain anyway.
# usage: waitchain.sh TAG    (TAG = z3u13-o18-lo on card 1, z3u13-o18-hi on card 2). Test hooks: ROOT, WPID, WSLEEP, DRY, NOPGREP.
cd ${ROOT:-/root/ecr} || exit 1; TAG=$1; O=out/z3u13-o18; E=enum/z3u13-o18; L=out/z3u12-o15r/run.log
log() { echo "$(date -u +%FT%TZ) $*" >> dside/chain.log; }
P=${WPID:-$(pgrep -x lotsieve7v | head -1)}
n0=$(wc -l < $L)
log "waitchain $TAG: waiting for lotsieve7v pid ${P:-none} (run.log has $n0 lines)"
while [ -n "$P" ] && kill -0 $P 2>/dev/null; do sleep ${WSLEEP:-10}; done
last=$(tail -n +$((n0 + 1)) $L | grep -E "^(DONE|STOPPED)" | tail -1 | cut -c1-80)
log "3/12 sweep exited; new status line: ${last:-none}"
case "$last" in
  STOPPED*) log "STOPPED: not chaining"; exit 0;;
  DONE*) ;;
  *) log "CRASH: sweep exited without a DONE/STOPPED line (run.log tail: $(tail -2 $L | tr '\n' ' ' | cut -c1-200)); chaining anyway so the GPU does not idle";;
esac
[ -e STOP ] && { log "STOP file present: not chaining"; exit 0; }
if [ -z "$NOPGREP" ] && { pgrep -x lotsieve7 >/dev/null || pgrep -x lotsieve7v >/dev/null; }; then log "another sieve is running: not chaining"; exit 0; fi
{ [ -s $E.info ] && [ -s $E.vrows64 ]; } || { log "enum $E missing: not chaining"; exit 1; }
[ -n "$DRY" ] && { log "DRY: would launch $TAG"; exit 0; }
mkdir -p $O
echo "$(date -u +%FT%TZ) build-dside chain ($TAG): 3/13 class exactly 18 distinct primes <= 97, any exponents, b <= 3.6e25 (omegalist128 97 3.6e25 18 minus 19 = 1,013,535 b, sha256 56ca9f96bab6a4af), this card's half $(cat dside/o18list.txt), H [17.87, 20.1997) (board 20.199607 + 1e-4), rule --vmul 2,3; enum: $(tail -1 $E.info)" >> $O/command.txt
setsid nohup taskset -c 0-22 dside/lotsieve7v --table search/lottery/tables/z3 --outdir $O --spans $E --cut1 10 --cut2 28 --cut3 40 --hshell 17.87 20.1997 \
  --tag $TAG --batch 2e8 --stop /root/ecr/STOP --family z3w --no-watchdog-i-know --vmul > $O/run.log 2>&1 < /dev/null &
sleep 20; pgrep -x lotsieve7v >/dev/null && log "launched $TAG: $(tail -1 $O/progress.tsv | cut -f1-10)" || log "launch FAILED: $(tail -3 $O/run.log)"
