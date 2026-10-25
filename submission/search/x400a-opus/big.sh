#!/bin/bash
# big.sh N thr1 thr2 tag
N=$1; t1=$2; t2=$3; tag=$4
cd "$(dirname "$0")"; mkdir -p o
for k in $(seq 0 13); do ./x4 $N $t1 $t2 99 14 $k > o/$tag.$k 2>o/$tag.$k.err & done; wait
cat o/$tag.? o/$tag.?? > o/$tag.all; rm -f o/$tag.? o/$tag.??
# drop already ranked
cat o/*.all.r 2>/dev/null | awk '{print $1" "$2}' | sort -u > o/ranked.keys
awk 'NR==FNR{K[$0]=1; next} !(($1" "$2) in K)' o/ranked.keys o/$tag.all | sort -k4,4nr > o/$tag.new
echo "$(date +%H:%M:%S) $tag sieved $(wc -l < o/$tag.all) new $(wc -l < o/$tag.new)" >> xlog.txt
./rankx.sh o/$tag.new 14 99
echo "$(date +%H:%M:%S) $tag ranked: $(awk '{print $5}' o/$tag.new.r | sort -n | uniq -c | tr '\n' ' ')" >> xlog.txt
