#!/bin/bash
# z3sel_apply.sh FILE [DROPFILE]  (owner, 20 Nov 2026)
# Applies the 3-isogeny Selmer bound (work/z3sel_filter.sh, search/selmer/z3sel.gp; 0 violations on 16,098 ranked curves)
# to a ranker queue file in place: a z3 line is kept only if its bound >= the lowest rank that scores at its height
# (boards of 20 Nov 07:20Z: 3/9 9.512795, 3/10 11.963673, 3/11 13.999595, 3/12 17.864161, 3/13 20.199607,
# 3/14 22.493691, 3/15 28.255216, 3/16 new). Non-z3 lines and lines whose bound errs are kept unchanged.
# Dropped lines are appended to DROPFILE (default FILE.selmer-dropped). Run only while the ranker owning FILE is stopped.
set -euo pipefail
R=.
F=${1:?file}; D=${2:-$F.selmer-dropped}
T=$(mktemp -d $R/work/z3sel/tmp/apply.XXXXXX); trap 'rm -rf "$T"' EXIT
# target per line from H (field 'H=x' in sieve lines, else field 7)
awk -v T="$T" '{
  h = -1; for (i = 1; i <= NF; i++) if ($i ~ /^H=/) { h = substr($i, 3) + 0; break }
  if (h < 0) h = $7 + 0
  if ($1 != "z3") { print > (T "/keep.other"); next }
  t = 16
  if (h < 28.255216) t = 15; if (h < 22.493691) t = 14; if (h < 20.199607) t = 13; if (h < 17.864161) t = 12
  if (h < 13.999595) t = 11; if (h < 11.963673) t = 10; if (h < 9.512795) t = 9
  print > (T "/grp." t)
}' "$F"
: > "$T/kept"
for g in "$T"/grp.*; do
  [ -e "$g" ] || continue
  t=${g##*.}
  $R/work/z3sel_filter.sh "$t" "$g" "$T/out.$t" 2>> "$T/stats"
  # strip the appended ' B3=<n>' (the filter moves the region tag to NF-1; put the line back as it was)
  sed -E 's/ B3=[^ ]+( R:[^ ]+)$/\1/; s/ B3=[^ ]+$//' "$T/out.$t" >> "$T/kept"
done
[ -e "$T/keep.other" ] && cat "$T/keep.other" >> "$T/kept"
# dropped = input lines not kept (match on the first three fields)
awk 'NR==FNR{k[$1" "$2" "$3]=1; next} !(($1" "$2" "$3) in k)' "$T/kept" "$F" >> "$D"
n0=$(wc -l < "$F"); n1=$(wc -l < "$T/kept")
cp "$T/kept" "$F"
echo "z3sel_apply $(basename $F): $n0 -> $n1 lines kept ($(cat $T/stats 2>/dev/null | tr '\n' ' ' | cut -c1-300))"
