#!/bin/bash
# lotrank.sh -- build-lottery ranking driver (own code; run as a systemd user unit from a frozen copy).
# Every 60 s: new lines of logs/lottery/*/survivors.txt with S(1e4) >= cut (rankcuts.txt) go to work/build-lottery/rank/pending.txt;
# the top 40 pending lines by S(1e4) go through lotrank.gp (1 core); outputs go to logs/lottery/<region>/ranked.txt and
# logs/lottery/ranked-all.txt.  r_lo >= 11, or r_lo == 10 with class-minimum H below the board: .ec via tools/mkec.py in
# found/lottery/, then ec_verify (novelty gate on) -> .verify.txt; certified rank >= 11 -> .CERTIFIED (read by lotnotify.sh);
# certified rank-10 undercut -> one AUTO line appended to helpers/build-lottery.status.
# Stop: touch work/build-lottery/STOP-RANK.
R=.; L=$R/search/lottery; W=$R/work/build-lottery/rank; GP=$R/repo/tools/pari/bin/gp
CORE=${CORE:-3}; CORES=${CORES:-$CORE}; mkdir -p $W $R/found/lottery; cd $W
declare -A TORS=([z6k]=6 [z5]=5 [leg]=4x2 [z3]=3) NEWC=([z6k]=11 [z5]=11 [leg]=11 [z3]=16)
# board heights to beat per (family, rank) (logs/board_parent_verify.json, 19 Nov: min(baseline, board_best)); rank >= 11 = new cell
declare -A BOARD=([z6k/8]=15.260014 [z6k/9]=18.445122 [z6k/10]=21.987758 [z5/8]=15.726734 [z5/9]=20.436404 [z5/10]=26.003803 [leg/8]=15.323158 [leg/9]=19.016457 [leg/10]=24.92637 [z3/8]=8.064684 [z3/9]=9.512795 [z3/10]=12.669669 [z3/11]=14.74411 [z3/12]=17.864161 [z3/13]=20.199607 [z3/14]=22.493691 [z3/15]=28.255216)
log() { echo "$(date -u +%FT%TZ) $*" >> $W/rank.log; }
log "lotrank start pid $$ cores $CORES"
# a batch interrupted by a restart goes back to pending
if [ -s $W/batch.in ] && [ ! -e $W/batch.done ]; then cat $W/batch.in >> $W/scored.txt; log "requeued $(wc -l < $W/batch.in) lines of an interrupted batch"; fi
while true; do
  [ -e $R/work/build-lottery/STOP-RANK ] && { log "STOP-RANK, exit"; exit 0; }
  date -u +%s > $W/heartbeat
  for d in $R/logs/lottery/*/; do
    id=$(basename $d); s=$d/survivors.txt; [ -s $s ] || continue
    off=$(cat $d/rank.offset 2>/dev/null || echo 0); n=$(wc -l < $s)
    if [ "$n" -gt "$off" ]; then
      sed -n "$((off+1)),${n}p" $s | awk -v id=$id -v cf=$L/rankcuts.txt 'BEGIN{while((getline l < cf)>0){split(l,a," "); if(a[1]!~/^#/) c[a[1]]=a[2]}
          cut = ""; for(k in c) if(index(id, k) == 1 && length(k) > length(best)) {best = k; cut = c[k]}}   # region-prefix cut wins over the family cut
        { cc = (cut != "") ? cut : (($1 in c) ? c[$1] : ""); if(cc != "" && $6+0 >= cc+0) print $0, "R:" id }' >> $W/pending.txt
      echo $n > $d/rank.offset
    fi
  done
  if [ ! -s $W/pending.txt ] && [ ! -s $W/scored.txt ] && [ -s $W/backlog.txt ]; then head -2000 $W/backlog.txt >> $W/pending.txt; tail -n +2001 $W/backlog.txt > $W/backlog.tmp; mv $W/backlog.tmp $W/backlog.txt; fi
  if [ -s $W/pending.txt ] || [ -s $W/scored.txt ]; then
    nw=$(echo $CORES | wc -w)
    if [ -s $W/pending.txt ]; then   # score: S(1e5) for up to 150 per worker of the unscored lines (highest S1e4 first), in parallel
      { awk '$4!="TESTCURVE_CLASS"{print $1, $2, $3}' $R/logs/lottery/ranked-all.txt; awk '{print $1, $2, $3}' $W/scored.txt 2>/dev/null; } > $W/seen.keys
      awk 'NR==FNR{seen[$1" "$2" "$3]=1; next} !(($1" "$2" "$3) in seen){seen[$1" "$2" "$3]=1; print}' $W/seen.keys $W/pending.txt 2>/dev/null | sort -k6,6gr > $W/pending.sorted
      NS=$((150 * nw)); head -$NS $W/pending.sorted > $W/score.in; tail -n +$((NS+1)) $W/pending.sorted > $W/pending.txt
      rm -f $W/score.part* $W/score.out.*; split -n r/$nw -d $W/score.in $W/score.part
      i=0; pids=; for c in $CORES; do pf=$W/score.part$(printf %02d $i); [ -s $pf ] || { i=$((i+1)); continue; }
        (echo "INFILE=\"$pf\"; OUTFILE=\"$W/score.out.$i\";"; cat $L/lotscore.frozen.gp) | taskset -c $c $GP -q -D parisizemax=500000000 >> $W/gp.err 2>&1 & pids="$pids $!"; i=$((i+1)); done
      for p in $pids; do wait $p; done
      cat $W/score.out.* >> $W/scored.txt 2>/dev/null
    fi
    # rank: the top 40 per worker by S(1e5)
    awk '{s=-99; for(i=1;i<=NF;i++) if($i ~ /^S5=/) s=substr($i,4)+0; print s "\t" $0}' $W/scored.txt | sort -t$'\t' -k1,1gr | cut -f2- > $W/scored.sorted
    NB=$((40 * nw)); head -$NB $W/scored.sorted > $W/batch.in; tail -n +$((NB+1)) $W/scored.sorted > $W/scored.txt
    [ -s $W/batch.in ] || { sleep 5; continue; }
    rm -f $W/batch.out $W/batch.done $W/batch.part* $W/batch.out.*
    split -n r/$nw -d $W/batch.in $W/batch.part
    i=0; pids=; for c in $CORES; do pf=$W/batch.part$(printf %02d $i); [ -s $pf ] || { i=$((i+1)); continue; }
      (echo "INFILE=\"$pf\"; OUTFILE=\"$W/batch.out.$i\"; T1=120; T3=900; TARGET=11;"; cat $L/lotrank.frozen.gp) > $W/batch.gp.$i
      taskset -c $c $GP -q -D parisizemax=1000000000 < $W/batch.gp.$i >> $W/gp.err 2>&1 & pids="$pids $!"; i=$((i+1)); done
    for p in $pids; do wait $p || log "gp worker rc $? (see gp.err)"; done
    cat $W/batch.out.* > $W/batch.out 2>/dev/null
    [ -s $W/batch.out ] || { log "empty batch output; batch kept in failed.txt"; cat $W/batch.in >> $W/failed.txt; sleep 60; continue; }
    cat $W/batch.out >> $R/logs/lottery/ranked-all.txt; touch $W/batch.done
    while read -r line; do
      set -- $line; fam=$1 p=$2 q=$3
      reg=$(echo "$line" | grep -o ' R:[^ ]*' | cut -c4-); [ -n "$reg" ] && [ -d $R/logs/lottery/$reg ] && echo "$line" >> $R/logs/lottery/$reg/ranked.txt
      [ "$4" = TESTCURVE_CLASS ] && { log "test-curve class dropped: $fam $p $q"; continue; }
      [ $# -lt 16 ] && continue
      lo=${9}; hfib=$7; hmin=${14}; [ "$lo" -ge 8 ] 2>/dev/null || continue
      # target cell: the highest rank k <= lo where this fibre scores (k >= 11: new cell; 8 <= k <= 10: fibre H below the board)
      nc=${NEWC[$fam]:-11}; k=0; for kk in $(seq $lo -1 8); do if [ $kk -ge $nc ]; then k=$kk; break; fi
        b=${BOARD[$fam/$kk]}; [ -n "$b" ] && awk -v a=$hfib -v b=$b 'BEGIN{exit !(a<b)}' && { k=$kk; break; }; done
      echo "$line" >> $R/logs/lottery/rank8plus.txt
      [ $k = 0 ] && { log "rank $lo above every board: $fam $p $q H $hfib (class min $hmin)"; continue; }
      cell=${TORS[$fam]}/$k; [ $k -ge $nc ] && cell=${TORS[$fam]}/$lo
      base=$R/found/lottery/${cell/\//-}-$fam-$p-$q; base=${base//-\-/-m}
      grep -F "$fam $p $q " $R/logs/lottery/ranked-all.txt | tail -1 > $W/one.txt
      python3 $R/tools/mkec.py $W/one.txt $fam $p $q $base.ec "build-lottery $fam p=$p q=$q region $reg cell $cell" >> $W/rank.log 2>&1 || { log "mkec failed $fam $p $q"; continue; }
      if [ $k -lt $nc ] && [ $k -lt $lo ]; then   # score in the lower cell k: keep k points
        awk -v k=$k '/^P /{n++; if(n>k) next} /^#CLAIM/{print "#CLAIM rank=" k; next} {print}' $base.ec > $base.ec.tmp && mv $base.ec.tmp $base.ec; fi
      (cd $R/repo && EC_GP_BIN=$R/repo/tools/pari/bin/gp timeout 1800 .venv-bench/bin/python -m ec_verify $base.ec) > $base.verify.txt 2>&1
      cr=$(grep -o 'certified rank >= [0-9]*' $base.verify.txt | grep -o '[0-9]*$' | head -1); hv=$(grep -o 'H = [0-9.]*' $base.verify.txt | head -1 | cut -c5-)
      vc=$(grep -o 'cell [0-9x]*/[0-9]*' $base.verify.txt | head -1 | cut -c6-)
      log "verify $base.ec: certified ${cr:-none} cell ${vc:-?} H ${hv:-?} (r_lo $lo, target $cell, class min H $hmin)"
      if [ -n "$cr" ] && [ "$cr" -ge $nc ] && [ $k -ge $nc ]; then
        echo "cell $vc H $hv classminH $hmin file $base.ec logdir $R/logs/lottery/$reg" > $base.CERTIFIED
      elif [ -n "$cr" ] && [ "$cr" -ge $k ] && [ $k -lt $nc ]; then
        b=${BOARD[$fam/$k]}
        if awk -v a=${hv:-99} -v b=$b 'BEGIN{exit !(a<b)}'; then
          echo "AUTO $(date -u +%FT%TZ) UNDERCUT ${vc:-$cell} H $hv board $b file $base.ec log $R/logs/lottery/$reg" >> $R/helpers/build-lottery.status
        else log "verified H $hv not below board $b for $cell: no AUTO line"; fi
      fi
    done < $W/batch.out
  else
    sleep 60
  fi
done
