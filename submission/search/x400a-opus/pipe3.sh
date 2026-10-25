#!/bin/bash
# pipe.sh N qlo qhi tag frac   (sieve rows q in [qlo,qhi], stage A S30000, rank top frac)
N=$1; qlo=$2; qhi=$3; tag=$4; frac=${5:-0.25}; thrS=${6:-2}
cd "$(dirname "$0")"; mkdir -p o
GP=${GP:-gp}
if [ ! -f o/$tag.sieve ]; then
for k in $(seq 0 13); do ./x4b $N 18 -50 99 14 $k $qlo $thrS 0.5 $qhi > o/$tag.$k 2>o/$tag.$k.err & done; wait
cat o/$tag.? o/$tag.?? > o/$tag.sieve; rm -f o/$tag.? o/$tag.??
fi
echo "$(date +%H:%M:%S) $tag sieved $(wc -l < o/$tag.sieve)" >> xlog.txt
cat o/*.r 2>/dev/null | awk '{print $1" "$2}' | sort -u > o/ranked.keys
awk 'NR==FNR{K[$0]=1; next} !(($1" "$2) in K) {print $1,$2,$7,$5}' o/ranked.keys o/$tag.sieve > o/$tag.A
for ((k=0;k<14;k++)); do
  awk -v k=$k 'NR%14==k {printf "%s[%s,%s,%s,%s]", (c++?",":"L=["), $1,$2,$3,$4} END{if(!c) printf "L=["; print "];"}' o/$tag.A > o/$tag.A$k.gp
  rm -f o/$tag.A$k.out; (echo "read(\"s3.gp\"); read(\"o/$tag.A$k.gp\"); s3(L,\"o/$tag.A$k.out\")" | $GP -q -s 1G > /dev/null 2>&1) &
done; wait
cat o/$tag.A*.out > o/$tag.s3; rm -f o/$tag.A*.gp o/$tag.A*.out
n=$(wc -l < o/$tag.s3); m=$(python3 -c "print(int($n*$frac))")
sort -k4,4nr o/$tag.s3 | head -$m > o/$tag.B
echo "$(date +%H:%M:%S) $tag stageA $n -> rank $m" >> xlog.txt
./rankx.sh o/$tag.B 14 99
echo "$(date +%H:%M:%S) $tag ranked: $(awk '{print $5}' o/$tag.B.r | sort -n | uniq -c | tr '\n' ' ')" >> xlog.txt
sort -k5,5nr -k3,3n o/$tag.B.r | head -5 >> xlog.txt
