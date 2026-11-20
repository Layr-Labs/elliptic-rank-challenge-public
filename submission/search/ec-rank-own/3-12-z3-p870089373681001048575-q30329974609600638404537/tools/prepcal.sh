#!/bin/bash
R=.; cd $R; Z=work/build-lottery/z4
for D in 1 -1 2 -2 6 -6 7 -7 15 -15 26 -26; do for arm in cls ctl; do Dn=${D/-/m}; id=z4c-$Dn-$arm; mkdir -p logs/lottery/$id
  [ -e work/build-lottery/enum/$id.info ] && continue
  echo "$(date -u +%FT%TZ) Z/4 calibration (owner note 5): family z4b, D = $D, arm $arm (cls = 47-smooth u <= 2.2e10/|D| with omega >= 5, coprime to D; ctl = as many uniform random u in [1, 2.2e10/|D|] coprime to D, seed 4000+D; work/build-lottery/z4/mkcal.py), s free, integral H [8, 12.6); cuts 5/18/30 from the frozen Z/4 rank 9-11 minima; gate: rank >= 9 yield per candidate cls vs ctl (GO >= 5x, STOP < 2x)" >> logs/lottery/$id/command.txt
  cmd="search/lottery/bin/lotenum work/build-lottery/enum/$id z4b 22000000000 8.0 12.6 -1 0 --qlist $Z/ucal_${D}_$arm.bin --dconst $D"
  echo "$(date -u +%FT%TZ) enum (pre): $cmd" >> logs/lottery/$id/command.txt
  bash -c "$cmd" 2>> logs/lottery/$id/enum.err
  echo "$id $(sed -n 3p work/build-lottery/enum/$id.info | cut -d' ' -f1-2)"
done; done; echo DONE
