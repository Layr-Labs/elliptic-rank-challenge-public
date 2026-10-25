#!/bin/bash
# usage: z4run3.sh Hlo Hmax Dmax outdir TLO THI
HL=$1; H=$2; DM=$3; OUT=$4; TL=$5; TH=$6
mkdir -p $OUT
python3 - $H $DM > $OUT/jobs.txt <<'PY'
import sys,math
H=float(sys.argv[1]); DM=int(sys.argv[2])
for D in range(1,DM+1):
    if any(D%(p*p)==0 for p in range(2,int(D**0.5)+1)): continue
    SM=int(3*math.exp(H)/D)+10
    W=max(50,int(400/D))
    for a in range(1,SM,W): print(D,a,min(SM,a+W))
PY
cat $OUT/jobs.txt | xargs -P 32 -L 1 sh -c "[ -f $OUT/c.\$0.\$1.done ] || { TLO=$TL THI=$TH HLO=$HL P3=10000 THR3=38 ./z4s_new $H 0 14 25 300 2000 \$0 \$1 \$2 > $OUT/c.\$0.\$1.txt 2>>$OUT/log && touch $OUT/c.\$0.\$1.done; }"
echo z4-finished $OUT
