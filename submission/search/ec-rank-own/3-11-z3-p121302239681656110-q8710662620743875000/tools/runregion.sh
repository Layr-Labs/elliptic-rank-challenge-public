#!/bin/bash
# runregion.sh ID FAMILY TABLEDIR HLO HHI CUT1 CUT2 CUT3 DAUTO [extra sieve args]
# One build-lottery region on the GPU: enum prefix work/build-lottery/enum/ID (.rows/.spans/.info), outdir logs/lottery/ID.
# Resumes from the checkpoint if one exists.  On exit 0 writes logs/lottery/ID/DONE and one logs/swept.tsv row.
# Exit codes are the sieve's (0 done, 10 stopped by STOP file/signal, 3 watchdog, 4 buffer, 5 memory).
R=.; L=$R/search/lottery
ID=$1 FAM=$2 TAB=$3 HLO=$4 HHI=$5 C1=$6 C2=$7 C3=$8 DA=$9; shift 9
CONIC=; case $FAM in *@*) CONIC=${FAM#*@}; FAM=${FAM%@*};; esac   # conic region: FAM = base@F1a,F1b,F1c,F2a,F2b,F2c (lattice forms)
PAIR=; case $FAM in *:pair:*) PAIR=${FAM#*:pair:}; FAM=${FAM%%:pair:*};; esac   # smooth-pair region: FAM = fam:pair:PLIST:OMIN[:REQ] (lotpair)
OUT=$R/logs/lottery/$ID; EN=$R/work/build-lottery/enum/$ID; mkdir -p $OUT; cd $R
BIN=$L/bin/lotsieve7; [ -n "$PAIR" ] && BIN=$L/bin/lotpair   # lotsieve6 dedups triples only (n = 3); lotsieve7 == lotsieve6 on n = 3 (identical survivors/dups, val/s6 s7)
RES=; [ -e $OUT/checkpoint ] && RES=--resume
CMD=(taskset -c 3 nice -n 10 $BIN --table $TAB --outdir $OUT --spans $EN --cut1 $C1 --cut2 $C2 --cut3 $C3 --hshell $HLO $HHI --duty-auto $DA --tag lot-$ID \
     --batch 5e7 --stop $R/work/build-lottery/STOP --code "bin:$(sha256sum $BIN | cut -c1-16)" $RES "$@")
CMD+=(--family $FAM); [ -n "$CONIC" ] && CMD+=(--conic $CONIC)
if [ -n "$PAIR" ]; then IFS=: read PL OM RQ <<< "$PAIR"; CMD+=(--plist $R/$PL --omin ${OM:-0}); [ -n "$RQ" ] && CMD+=(--require $RQ); fi
echo "$(date -u +%FT%TZ) region $ID: $(head -1 $EN.info | tr '\n' ' ') $(sed -n 3p $EN.info)" >> $OUT/command.txt
echo "$(date -u +%FT%TZ) ${CMD[*]}" >> $OUT/command.txt
"${CMD[@]}" 2>> $OUT/stderr.txt
rc=$?
echo "$(date -u +%FT%TZ) exit $rc" >> $OUT/command.txt
if [ $rc = 0 ]; then
  date -u +%FT%TZ > $OUT/DONE
  cand=$(grep -o 'cand_run=[0-9]*' $OUT/params.txt | cut -d= -f2 | paste -sd+ | bc); nfin=$(wc -l < $OUT/survivors.txt)
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$(date -u +%FT%TZ)" "$FAM" "lottery:$ID:H[$HLO,$HHI)" "$(sed -n 1p $EN.info | tr ' ' ':')" \
    "lotsieve:S300>=$C1,S2000>=$C2,S10000>=$C3" "$cand" "$nfin" "bin:$(sha256sum $BIN | cut -c1-16)" >> $R/logs/swept.tsv
fi
exit $rc
