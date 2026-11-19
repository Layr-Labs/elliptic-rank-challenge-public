#!/usr/bin/env python3
"""rowsfilter.py PREFIX MOD -- keep only rows whose outer value is coprime to MOD (z3 with 2^e 3^f | b: rows with 2 | a or 3 | a are
never weighted-primitive), rewriting PREFIX.rows/.spans in place (via temp files + rename) and appending the new counts to PREFIX.info."""
import sys, struct, math, os
pfx, mod = sys.argv[1], int(sys.argv[2])
r64 = os.path.exists(pfx + '.rows64'); ext = '.rows64' if r64 else '.rows'; rs = 24 if r64 else 16
rb = open(pfx + ext, 'rb').read(); sb = open(pfx + '.spans', 'rb').read()
nr = len(rb) // rs; keep = {}; out_r = bytearray()
for k in range(nr):
    q = struct.unpack_from('<q', rb, rs * k)[0]
    if math.gcd(q, mod) == 1: keep[k] = len(keep); out_r += rb[rs * k:rs * k + rs]
out_s = bytearray(); pos = 0
for k in range(len(sb) // 12):
    r, i0, ln = struct.unpack_from('<iii', sb, 12 * k)
    if r in keep: out_s += struct.pack('<iii', keep[r], i0, ln); pos += ln
open(pfx + ext + '.tmp', 'wb').write(out_r); open(pfx + '.spans.tmp', 'wb').write(out_s)
os.rename(pfx + ext + '.tmp', pfx + ext); os.rename(pfx + '.spans.tmp', pfx + '.spans')
open(pfx + '.info', 'a').write('rowsfilter gcd(outer,%d)=1: rows %d -> %d, spans %d, positions %.6e\n' % (mod, nr, len(keep), len(out_s) // 12, pos))
print('rows', nr, '->', len(keep), 'positions', '%.4g' % pos)
