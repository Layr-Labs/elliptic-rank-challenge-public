#!/bin/bash
# watch_ec.sh - bounded (<= 28 min) watcher run from the Mac: exits early when a curve with rank >= 9 appears anywhere in the Spark campaigns
# (found/*.ec, FOUND.txt, a ranked line with lb >= 9 or ub >= 9); otherwise prints a status digest after 28 min.  Re-armed by hand, never loops forever.
end=$((SECONDS+1680))
SSH="ssh -o BatchMode=yes -o ConnectTimeout=20 -J myspark2 spark-new"
while [ $SECONDS -lt $end ]; do
  out=$($SSH 'cd ~/ecwork_src/kub; for d in runs/camp_*; do ls $d/found/*.ec 2>/dev/null | sed "s/^/EC /"; [ -f $d/FOUND.txt ] && echo "FOUNDTXT $d"; done; cat runs/camp_*/inc/ranked_*.txt runs/camp_*/ranked.txt 2>/dev/null | awk "NF>=10 && \$7 ~ /^-?[0-9]+\$/ && (\$7>=9 || \$8>=9)" | sort -u | sed "s/^/HI /" | cut -c1-200; cat runs/camp_88*/inc/ranked_*.txt runs/camp_88*/ranked.txt 2>/dev/null | awk "NF>=10 && \$7 ~ /^-?[0-9]+\$/ && (\$7>=7 || \$8>=8)" | sort -u | sed "s/^/Z8 /" | cut -c1-200' 2>/dev/null)
  if [ -n "$out" ]; then echo "$out" | head -20; echo "WATCH: candidate event at $(date -u +%T)Z"; exit 0; fi
  sleep 60
done
$SSH 'cd ~/ecwork_src/kub; date -u +%T; for f in campaign_f44qa3.log campaign_f44qb3.log campaign_z5qa3.log campaign_z5qb3.log campaign_z66q3.log campaign_z67q3.log campaign_z88q3.log campaign_f44qa2p.log campaign_f44qb2p.log campaign_z5qa2p.log campaign_z5qb2p.log campaign_z5qa4.log campaign_z5qb4.log; do echo "-- $f"; grep "== shell" $f | tail -2 | cut -c1-120; done; ./hi_scan.sh 44 | tail -8 | cut -c1-170; ./hi_scan.sh 6 | tail -8 | cut -c1-170; ps -eo pid,stat,args | grep "kubg" | grep -v grep | cut -c1-100' 2>&1
echo "WATCH: quiet (no find) at $(date -u +%T)Z"
