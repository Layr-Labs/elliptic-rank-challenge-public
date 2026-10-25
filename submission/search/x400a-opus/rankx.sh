#!/bin/bash
# rankx.sh infile nshards Hmax -> infile.r
in=$1; ns=${2:-14}; H=${3:-30}
GP=${GP:-gp}; D=$(cd "$(dirname "$0")" && pwd)
for ((k=0;k<ns;k++)); do
  awk -v k=$k -v n=$ns 'NR%n==k' $in | awk 'BEGIN{printf "L=["} {printf "%s[%s,%s]", (NR>1?",":""), $1,$2} END{print "];"}' > $in.s$k.gp
  rm -f $in.s$k.r
  (echo "read(\"$D/lib.gp\"); read(\"$in.s$k.gp\"); st(L,$H,0,\"$in.s$k.r\");" | $GP -q -s 2G > /dev/null 2>&1) &
done; wait
cat $in.s*.r > $in.r; rm -f $in.s*.gp $in.s*.r
