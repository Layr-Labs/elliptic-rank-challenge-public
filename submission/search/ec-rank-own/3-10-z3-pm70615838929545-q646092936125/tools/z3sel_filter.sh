#!/bin/bash
# z3sel_filter.sh -- drop Z/3 Tate-family candidates whose 3-isogeny Selmer bound is below a target rank.
#   usage: work/z3sel_filter.sh TARGET IN [OUT]          (OUT defaults to stdout; stats go to stderr)
#   IN: sieve/pending lines 'z3 P Q ...' with P = b, Q = a^3 (curve y^2 + a x y + b y = x^3); any other first field
#   is passed through untouched with ' B3=na' (not this family). Kept z3 lines get ' B3=<bound>' appended, where
#   bound = s_phi + s_phihat - 1 >= rank (search/selmer/z3sel.gp, closed form, validated against oracle3.gp).
#   A line whose bound cannot be computed is KEPT with ' B3=err' (fail-safe).  The region tag moves from $NF to $(NF-1).
#   One gp process, pinned to core 19, capped at 2 GB.  Work files under work/z3sel/tmp.
set -euo pipefail
ROOT=.
[ $# -ge 2 ] || { echo "usage: $0 TARGET IN [OUT]" >&2; exit 2; }
TGT=$1; IN=$2; OUT=${3:-/dev/stdout}
[[ "$TGT" =~ ^[0-9]+$ ]] || { echo "TARGET must be a non-negative integer" >&2; exit 2; }
[ -r "$IN" ] || { echo "cannot read $IN" >&2; exit 2; }
mkdir -p "$ROOT/work/z3sel/tmp"
TMP=$(mktemp -d "$ROOT/work/z3sel/tmp/f.XXXXXX")
trap 'rm -rf "$TMP"' EXIT
t0=$(date +%s.%N)
cp "$IN" "$TMP/in.txt"       # read the input once: the live rank files are rewritten while the ranker runs
# z3 lines -> [lineno, P, Q]
awk '$1=="z3" && $2 ~ /^-?[0-9]+$/ && $3 ~ /^-?[0-9]+$/ {print "[" NR "," $2 "," $3 "]"}' "$TMP/in.txt" > "$TMP/in.gpv"
cat > "$TMP/run.gp" <<EOF
read("$ROOT/search/selmer/z3sel.gp");
V = readvec("$TMP/in.gpv");
{ for(i = 1, #V, my(v = V[i], r = iferr(z3bound(v[2], v[3]), e, 0));
    print(v[1], " ", v[2], " ", v[3], " ", if(r == 0, "err", r[1]))); }
quit
EOF
systemd-run --user --scope -q -p MemoryMax=2G taskset -c 19 nice -n 15 "$ROOT/repo/tools/pari/bin/gp" -q -f "$TMP/run.gp" > "$TMP/out.txt" < /dev/null
nz=$(wc -l < "$TMP/in.gpv"); nb=$(wc -l < "$TMP/out.txt")
[ "$nz" -eq "$nb" ] || { echo "z3sel_filter: gp returned $nb bounds for $nz z3 lines -- aborting, nothing written" >&2; exit 1; }
awk -v T="$TGT" -v S="$TMP/stats" 'NR==FNR {B[$1]=$4; KP[$1]=$2; KQ[$1]=$3; next}
  { if (!(FNR in B)) { if ($1=="z3") { print $0 " B3=err"; ne++ } else { print $0 " B3=na"; nn++ }; next }
    if (($2 "") != (KP[FNR] "") || ($3 "") != (KQ[FNR] "")) { print "z3sel_filter: join mismatch at line " FNR > "/dev/stderr"; bad=1; exit 1 }
    nzl++;
    if (B[FNR]=="err") { print $0 " B3=err"; ne++; next }
    if (B[FNR]+0 >= T) { print $0 " B3=" B[FNR]; nk++ } else nd++ }
  END { printf "z3 lines %d  kept %d  dropped %d  errors(kept) %d  non-z3 passed %d\n", nzl, nk, nd, ne, nn > S }' \
  "$TMP/out.txt" "$TMP/in.txt" > "$TMP/kept.txt"
cat "$TMP/kept.txt" > "$OUT"
t1=$(date +%s.%N)
echo "z3sel_filter target=$TGT in=$IN: $(cat "$TMP/stats")  ($(echo "$t1 - $t0" | bc) s)" >&2
