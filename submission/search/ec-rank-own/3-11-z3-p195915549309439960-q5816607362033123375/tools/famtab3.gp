\\ famtab3.gp -- own a_l tables for the build-lottery families (generalises search/gpusieve/z6/z6tab.gp).
\\ Needs fam3.gp (famint).  Index t = P * Q^-1 mod l (0 <= t < l), index l = class Q == 0 mod l.
\\ Every a_i(P,Q) of famint is homogeneous of degree i in (P,Q), so for l not dividing Q the reduction is the model at
\\ (t,1) scaled by u = Q: a_l depends on t only.  Classes with Delta(t) == 0 mod l: if c4 != 0 mod l the reduction is
\\ multiplicative and the model is minimal at l, a_l = kronecker(-c6, l); if c4 == 0 too (additive), a_l = 0 and the
\\ class is flagged ADDITIVE (checked by lifts like every class).
\\ Checks: every class for l < LIFTALL, else 3 random classes + all bad classes, each against ellap of the GLOBAL integral
\\ curve at 2 random coprime lifts (|P|, |Q| ~ 1e12, mixed signs; leg: positive).
LIFTALL = 500;
cvec(f, P, Q) = { my(a = famint(f, P, Q), a1 = a[1], a2 = a[2], a3 = a[3], a4 = a[4], a6 = a[5],
  b2 = a1^2+4*a2, b4 = 2*a4+a1*a3, b6 = a3^2+4*a6, b8 = a1^2*a6+4*a2*a6-a1*a3*a4+a2*a3^2-a4^2,
  c4 = b2^2-24*b4, c6 = -b2^3+36*b2*b4-216*b6, D = -b2^2*b8-8*b4^3-27*b6^2+9*b2*b4*b6); [c4, c6, D] };
okpair(f, P, Q) = gcd(P,Q) == 1 && P != 0 && Q != 0 && cvec(f, P, Q)[3] != 0 && (f != "leg" || (P > 0 && Q > 0));
liftpq(f, t, l) = { my(P, Q, sg = (f != "leg"));
  until(okpair(f, P, Q),
    if(t < l, Q = (random(10^12)+1)*if(sg, 2*random(2)-1, 1); while(Q % l == 0, Q++); P = lift(Mod(t,l)*Q) + l*(random(2*10^12)-10^12); if(!sg, P = abs(P)); if(!sg && P % l != lift(Mod(t,l)*Q), P = 0),
              Q = l*(random(10^11)+1)*if(sg, 2*random(2)-1, 1); P = (random(10^12)+1)*if(sg, 2*random(2)-1, 1)));
  [P,Q] };
run(f, lo, hi, out, lg) = {
  my(nbad = 0, nadd = 0, nchk = 0, nmis = 0, t0 = getwalltime());
  forprime(l = max(lo,5), hi-1,
    my(v = vector(l+1), bad = List());
    for(t = 0, l,
      my(pq = if(t < l, [t,1], [1,0]), cv = cvec(f, pq[1], pq[2]), a);
      if(cv[3] % l != 0, a = ellap(ellinit(famint(f, pq[1], pq[2]), l)),
         listput(bad, t); nbad++;
         if(cv[1] % l != 0, a = kronecker(-cv[2], l), a = 0; nadd++; write(lg, "ADDITIVE l=", l, " t=", t)));
      v[t+1] = a);
    my(chk = if(l < LIFTALL, [0..l], concat(vector(3, j, random(l+1)), Vec(bad))));
    if(f == "z3", chk = select(t -> t < l, chk));   \\ z3: class Q = 0 (l | a) is 0 by convention (true a_l depends on b mod l), not lift-checked
    for(j = 1, #chk, my(t = chk[j], pq = liftpq(f, t, l), b = ellap(famE(f, pq[1], pq[2]), l), pq2 = liftpq(f, t, l), b2 = ellap(famE(f, pq2[1], pq2[2]), l));
      nchk++; if(b != v[t+1] || b2 != v[t+1], nmis++; write(lg, "MISMATCH l=", l, " t=", t, " table ", v[t+1], " lifts ", [b,b2], " pq ", [pq, pq2])));
    write(out, Str(l, " ", strjoin(apply(x->Str(x), v), " "))));
  write(lg, f, " range [", lo, ",", hi, ") bad classes ", nbad, " additive ", nadd, " lift checks ", nchk, " mismatches ", nmis, " ms ", getwalltime()-t0);
};
