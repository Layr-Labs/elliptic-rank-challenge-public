#!/bin/bash
# hi_scan.sh - list every ranked curve of the campaigns with lower bound >= 8 or upper bound >= 9 (candidates for a deeper point search), plus any found/*.ec and FOUND.txt.
# Runs on the Spark from ~/ecwork_src/kub:  ./hi_scan.sh [family-substring, default all]
cd ~/ecwork_src/kub
pat=${1:-}
for d in runs/camp_${pat}*; do
  [ -d "$d" ] || continue
  f=$(cat "$d"/inc/ranked_*.txt "$d"/ranked.txt 2>/dev/null | awk 'NF>=10 && $7 ~ /^-?[0-9]+$/ && ($7>=8 || $8>=9)' | sort -u)
  if [ -n "$f" ]; then echo "== $d"; echo "$f" | cut -c1-160; fi
  ls "$d"/found/*.ec 2>/dev/null | sed 's/^/EC FILE: /'
  [ -f "$d/FOUND.txt" ] && { echo "FOUND.txt in $d:"; cut -c1-250 "$d/FOUND.txt"; }
done
