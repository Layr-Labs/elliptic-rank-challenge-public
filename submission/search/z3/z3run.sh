#!/bin/bash
# usage: z3run.sh Hlo Hmax outdir [chunkA]   (chunked, resumable)
HL=$1; H=$2; OUT=$3; CA=${4:-10}
mkdir -p $OUT
python3 -c "
import math
AM=int(math.exp($H+0.06))+2
for a in range(1,AM,$CA): print(a,min(AM,a+$CA))" > $OUT/jobs.txt
cat $OUT/jobs.txt | xargs -P 32 -L 1 sh -c "[ -f $OUT/c.\$0.done ] || { HLO=$HL P3=10000 THR3=${THR3:-38} ./z3s_new $H 0 ${THR1:-15} ${THR2:-28} 300 2000 \$0 \$1 > $OUT/c.\$0.txt 2>>$OUT/log && touch $OUT/c.\$0.done; }"
echo z3-finished $OUT
