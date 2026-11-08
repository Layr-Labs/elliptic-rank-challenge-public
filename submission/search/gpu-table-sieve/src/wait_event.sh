#!/bin/bash
# wait_event.sh [seconds=570] - foreground wait (<= the given time) that returns early when a curve with rank lower bound >= 9 or upper bound >= 9 shows up in any campaign (or a found/*.ec appears).
# Prints the event lines, or a one-line digest when the time is up.
end=$((SECONDS+${1:-570}))
SSH="ssh -o BatchMode=yes -o ConnectTimeout=20 -J myspark2 spark-new"
while [ $SECONDS -lt $end ]; do
  out=$($SSH 'cd ~/ecwork_src/kub; for d in runs/camp_*; do ls $d/found/*.ec 2>/dev/null | sed "s/^/EC /"; [ -f $d/FOUND.txt ] && echo "FOUNDTXT $d"; done; cat runs/camp_*/inc/ranked_*.txt runs/camp_*/ranked.txt 2>/dev/null | awk "NF>=10 && \$7 ~ /^-?[0-9]+\$/ && (\$7>=9 || \$8>=9)" | sort -u | sed "s/^/HI /" | cut -c1-200' 2>/dev/null)
  if [ -n "$out" ]; then echo "$out" | head -20; echo "EVENT at $(date -u +%T)Z"; exit 0; fi
  sleep 45
done
$SSH 'cd ~/ecwork_src/kub; echo "quiet at $(date -u +%T)Z; running sieves: $(ps -eo args | grep -c "[s]ieve [0-9]")"; for f in campaign_f44qa3 campaign_z5qa3 campaign_z66q3 campaign_z67q3 campaign_z5qa4 campaign_f44qa4; do echo "$f: $(grep "== shell" $f.log | tail -1 | cut -c1-60)"; done' 2>&1
