\\ fam3.gp -- the three families of build-lottery (own code).  Integral models in the sieve coordinates:
\\   z6k: Kubert Z/6 in (pk,qk) = (p, p+q), c = pk/(qk-pk):  [qk-2pk, -pk qk, -pk qk (qk-pk), 0, 0]
\\   z5 : Kubert Z/5 (Tate normal form b = c = t, t = p/q):  [q-p, -p q, -p q^2, 0, 0]
\\   leg: Legendre Z/2xZ/4, y^2 = x(x+a^2)(x+b^2) in (p,q) = (a,b):  [0, a^2+b^2, 0, a^2 b^2, 0]
\\   z3 : Z/3, y^2 + a x y + b y = x^3 ~ E_s: y^2 + x y + s y = x^3 with s = b/a^3 = P/Q:  [Q, 0, P Q^2, 0, 0]
\\        (famcurve(t) = [1, 0, t, 0, 0] with t = p/q; the sieve writes p = b, q = a^3, i.e. t = b/a^3)
\\   z4 : Z/4, Tate y^2 + x y - t y = x^3 - t x^2, t = P/Q (the sieve writes P = u, Q = D s^2):  [Q, -P Q, -P Q^2, 0, 0]
\\ S(L) = sum_{5 <= l < L} -a_l log(l)/l (the sieve convention: primes 2, 3 left out).
famint(f, P, Q) = {
  if(f == "z6k", [Q-2*P, -P*Q, -P*Q*(Q-P), 0, 0],
     f == "z5",  [Q-P, -P*Q, -P*Q^2, 0, 0],
     f == "leg", [0, P^2+Q^2, 0, P^2*Q^2, 0],
     f == "z3",  [Q, 0, P*Q^2, 0, 0],
     f == "z4",  [Q, -P*Q, -P*Q^2, 0, 0], error("family ", f)) };
famE(f, P, Q) = ellinit(famint(f, P, Q));
heightE(E) = log(max(abs(E.c4)^3, E.c6^2))/12;
famH(f, P, Q) = heightE(famE(f, P, Q));
famHmin(f, P, Q) = heightE(ellminimalmodel(famE(f, P, Q)));
\\ sums at the sieve's stage limits (from 5) for an ellinit'd curve, vector over lims
Ssums(E, lims) = { my(a = ellan(E, lims[#lims]), S = 0., r = vector(#lims), k = 1);
  forprime(l = 5, lims[#lims], while(l >= lims[k], r[k] = S; k++); S -= a[l]*log(l)/l); r[#lims] = S; r };
Wp = [7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47];
wcount(n) = sum(j = 1, #Wp, n % Wp[j] == 0);
wset(n) = select(w -> n % w == 0, Wp);
