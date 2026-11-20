\\ z6tab.gp -- own a_l table for the Kubert Z/6 family in (pk,qk) = (p, p+q) coordinates (c = pk/(qk-pk)).
\\ Integral model E(pk,qk) = [qk-2pk, -pk qk, -pk qk (qk-pk), 0, 0], Delta = pk^6 qk^3 (qk+8pk) (qk-pk)^2.
\\ Index t = pk * qk^-1 mod l (0 <= t < l), index l = class qk == 0 mod l.
\\ Good classes: a_l of the reduction of [1-2t, -t, -t(1-t), 0, 0] mod l (t < l); a_l(pk,qk) depends only on t since
\\   a_i(pk,qk) = qk^i a_i(t,1) mod l (weighted homogeneous, u = qk is a unit).
\\ Bad classes t = 0, 1, -1/8, infinity: c4 = 1, 9, 9/256, 16 there (nonzero mod l >= 5), so the reduction is
\\   multiplicative and the model is minimal at l: a_l = kronecker(-c6, l).
\\ Checks: every class for l < LIFTALL and 3 random classes per larger prime (all bad classes always) against
\\   ellap of the GLOBAL integral curve at 2 random coprime lifts (|pk|, |qk| ~ 1e12, mixed signs).
LIFTALL = 500;
E6(P,Q) = ellinit([Q-2*P, -P*Q, -P*Q*(Q-P), 0, 0]);
c6of(P,Q) = { my(a1 = Q-2*P, a2 = -P*Q, a3 = -P*Q*(Q-P), b2 = a1^2+4*a2, b4 = a1*a3, b6 = a3^2); -b2^3 + 36*b2*b4 - 216*b6 };
liftpq(t, l) = { my(P, Q);
  until(gcd(P,Q) == 1 && P != 0 && Q != 0 && Q != P && Q != -8*P,
    if(t < l, Q = (random(10^12)+1)*(2*random(2)-1); while(Q % l == 0, Q++); P = lift(Mod(t,l)*Q) + l*(random(2*10^12)-10^12),
              Q = l*(random(10^11)+1)*(2*random(2)-1); P = (random(10^12)+1)*(2*random(2)-1)));
  [P,Q] };
run(lo, hi, out, lg) = {
  my(nbad = 0, nchk = 0, nmis = 0, t0 = getwalltime());
  forprime(l = max(lo,5), hi-1,
    my(v = vector(l+1), chk = if(l < LIFTALL, [0..l], concat(vector(3, j, random(l+1)), [0, 1, lift(Mod(-1,l)/8), l])));
    for(t = 0, l,
      my(a);
      if(t < l && (t * (1+8*t) * (1-t)) % l != 0,
        a = ellap(ellinit([1-2*t, -t, -t*(1-t), 0, 0], l)),
        my(pq = if(t < l, [t,1], [1,0])); a = kronecker(-c6of(pq[1], pq[2]), l); nbad++);
      v[t+1] = a);
    for(j = 1, #chk, my(t = chk[j], pq = liftpq(t, l), b = ellap(E6(pq[1], pq[2]), l), pq2 = liftpq(t, l), b2 = ellap(E6(pq2[1], pq2[2]), l));
      nchk++; if(b != v[t+1] || b2 != v[t+1], nmis++; write(lg, "MISMATCH l=", l, " t=", t, " table ", v[t+1], " lifts ", [b,b2])));
    write(out, Str(l, " ", strjoin(apply(x->Str(x), v), " "))));
  write(lg, "range [", lo, ",", hi, ") bad classes ", nbad, " lift checks ", nchk, " mismatches ", nmis, " ms ", getwalltime()-t0);
};
