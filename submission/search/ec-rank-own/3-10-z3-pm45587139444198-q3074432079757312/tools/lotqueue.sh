#!/bin/bash
# lotqueue.sh -- build-lottery GPU queue driver (own code; systemd user unit lottery-queue, frozen copy).
# Re-reads search/lottery/queue.tsv before every region (lines: id fam table hlo hhi cut1 cut2 cut3 dauto enum-command...;
# '#' = comment); runs the first region without logs/lottery/<id>/DONE: enumerates it if work/build-lottery/enum/<id>.info is
# missing, then bin/runregion.frozen.sh (resumes from its checkpoint).  rc 0 -> next region (enum files removed, .info kept);
# rc 10 -> STOP file present: exit 0, else retry after 60 s; other rc -> logs/lottery/<id>/FAIL, exit with that rc.
# Queue exhausted -> work/build-lottery/QUEUE-EMPTY, exit 0.  Stop cleanly: touch work/build-lottery/STOP (the sieve
# checkpoints and exits 10).
R=.; L=$R/search/lottery; WB=$R/work/build-lottery; LOG=$WB/queue.log
cd $R; rm -f $WB/QUEUE-EMPTY
log() { echo "$(date -u +%FT%TZ) $*" >> $LOG; }
log "lotqueue start pid $$"
while systemctl --user is-active --quiet lottery-gpu-z6-r1; do date -u +%s > $WB/queue.heartbeat; sleep 30; done   # (first-region unit of 07:23Z; gone after the 07:47Z reorder)
while true; do
  [ -e $WB/STOP ] && { log "STOP present, exit"; exit 0; }
  date -u +%s > $WB/queue.heartbeat
  next=$(grep -v '^#' $L/queue.tsv | awk 'NF>=10' | while read -r id rest; do [ -e $R/logs/lottery/$id/DONE ] || { echo "$id $rest"; break; }; done)
  [ -z "$next" ] && { log "queue empty"; date -u +%FT%TZ > $WB/QUEUE-EMPTY; exit 0; }
  set -- $next; id=$1 fam=$2 tab=$3 hlo=$4 hhi=$5 c1=$6 c2=$7 c3=$8 da=$9; shift 9; ecmd="$*"
  if [ ! -e $WB/enum/$id.info ]; then
    log "enum $id: $ecmd"; mkdir -p $R/logs/lottery/$id; echo "$(date -u +%FT%TZ) enum: $ecmd" >> $R/logs/lottery/$id/command.txt
    taskset -c 3 nice -n 10 bash -c "$ecmd" >> $R/logs/lottery/$id/enum.err 2>&1
    [ -e $WB/enum/$id.info ] || { log "enum failed for $id"; echo "enum failed" > $R/logs/lottery/$id/FAIL; exit 6; }
  fi
  log "run $id"
  $L/bin/runregion.frozen.sh $id $fam $tab $hlo $hhi $c1 $c2 $c3 $da
  rc=$?; log "region $id rc $rc"
  case $rc in
    0) rm -f $WB/enum/$id.rows $WB/enum/$id.spans;;
    10) [ -e $WB/STOP ] && { log "stopped by STOP, exit"; exit 0; }; sleep 60;;
    5) log "region $id: not enough free memory (other lanes); retry in 10 min"; sleep 600;;
    *) echo "rc $rc $(date -u +%FT%TZ)" > $R/logs/lottery/$id/FAIL; exit $rc;;
  esac
done
