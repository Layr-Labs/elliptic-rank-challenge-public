#!/usr/bin/env python3
"""build-lottery tables: part*.txt ('l a_0 .. a_l') -> al.bin (int16 LE, l+1 entries per prime) + al.json, the schema of
search/gpusieve/tables/z6k (index offsets[i] + t: class t = P Q^-1 mod l; offsets[i] + l: class Q = 0 mod l).
usage: al2bin3.py DIR FAMILY 'coords text' part1.txt part2.txt ..."""
import sys, json, array, hashlib, os
d, fam, coords = sys.argv[1:4]; rows = {}
for f in sys.argv[4:]:
    for line in open(os.path.join(d, f)):
        v = line.split()
        if v: l = int(v[0]); a = [int(x) for x in v[1:]]; assert len(a) == l + 1 and all(x * x <= 4 * l for x in a), l; rows[l] = a
isprime = lambda n: n > 1 and all(n % k for k in range(2, int(n ** .5) + 1))
need = [p for p in range(5, 10000) if isprime(p)]
missing = [p for p in need if p not in rows]; assert not missing, missing[:5]
buf = array.array('h'); offs = []
for l in need: offs.append(len(buf)); buf.extend(rows[l])
assert sys.byteorder == 'little'
open(os.path.join(d, 'al.bin'), 'wb').write(buf.tobytes())
json.dump({'family': fam, 'coords': coords, 'index': 'offsets[i] + t: class t = P * Q^-1 mod l; offsets[i] + l: class Q = 0 mod l',
           'primes': need, 'offsets': offs, 'complete_upto': 10000, 'dtype': 'int16', 'al_bin_sha256': hashlib.sha256(buf.tobytes()).hexdigest()},
          open(os.path.join(d, 'al.json'), 'w'))
print(fam, 'primes', len(need), 'entries', len(buf))
