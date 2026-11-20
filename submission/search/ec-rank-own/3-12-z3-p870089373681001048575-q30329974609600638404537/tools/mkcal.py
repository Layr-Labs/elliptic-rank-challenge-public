#!/usr/bin/env python3
"""Z/4 calibration lists (build-lottery note 5): per D, class = 47-smooth u <= 2.2e10/|D| with omega >= 5 (u47k5.bin), coprime to D;
control = the same number of uniform random u in [1, 2.2e10/|D|], coprime to D (seed 4000 + D).  Writes ucal_<D>_{cls,ctl}.bin (int64, sorted)."""
import array, random, math
U = array.array('q'); U.frombytes(open('u47k5.bin', 'rb').read())
for D in (1, -1, 2, -2, 6, -6, 7, -7, 15, -15, 26, -26):
    X = 22000000000 // abs(D)
    cls = [u for u in U if u <= X and math.gcd(u, abs(D)) == 1]
    rnd = random.Random(4000 + D); ctl = set()
    while len(ctl) < len(cls):
        u = rnd.randint(1, X)
        if math.gcd(u, abs(D)) == 1: ctl.add(u)
    open('ucal_%d_cls.bin' % D, 'wb').write(array.array('q', sorted(cls)).tobytes())
    open('ucal_%d_ctl.bin' % D, 'wb').write(array.array('q', sorted(ctl)).tobytes())
    print(D, len(cls))
