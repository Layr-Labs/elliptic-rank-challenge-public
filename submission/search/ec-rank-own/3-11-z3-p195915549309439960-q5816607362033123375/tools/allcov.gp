\\ allcov.gp -- 2-coverings of EVERY class of Sel2(E)/delta(E(Q)_tors), for any E over Q (proof-covers brief).
\\ Own code (GP). Model y^2 = f(x) = x^3 + b2 x^2 + 8 b4 x + 16 b6 (change v = [1/2, 0, -a1/2, -a3/2] from E).
\\ L = Q[th]/(f) (any factorisation of f). Selmer element delta in L*/L*^2:
\\   from an ell2cover cover [R, [X, Y]] (X = N(t)/y^2, y^2 = R(t)):  4 N(t) - th R(t) = delta * (quadratic in t)^2 in L[t];
\\   from a point P: x(P) - th (componentwise fix at a 2-torsion point).
\\ Cover of a class delta: delta (u0 + u1 th + u2 th^2)^2 = Q0(u) + Q1(u) th + Q2(u) th^2; conic Q2(u) = 0 parametrised
\\ u = Mp [s^2, s t, t^2]~; quartic w^2 = -Q1(u(s,t)); the point is x = -Q0(u)/Q1(u) on y^2 = f(x).
\\ F2 linear algebra on classes by quadratic characters at degree-1 primes (l, r) of L (r root of f mod l).

fmodel(E) = { my(v = [1/2, 0, -E.a1/2, -E.a3/2], f = 'x^3 + E.b2*'x^2 + 8*E.b4*'x + 16*E.b6, E2 = ellchangecurve(E, v));
  if(E2[1..5] != [0, E.b2, 0, 8*E.b4, 16*E.b6], error("fmodel"));
  [v, f] };
\\ L-element from a polynomial in th
Lel(f, p) = Mod(subst(p, 'x, 'th) * 1, subst(f, 'x, 'th));
Lnorm(d) = { my(f = d.mod, a = lift(d)); polresultant(f, a, 'th) };   \\ f monic

\\ delta of an ell2cover cover (on E as given to ell2cover); f from fmodel(E)
\\ c * square test for a polynomial over a field (polmod coefficients allowed), degree 0, 2 or 4
iscsquare(P) = {
  my(d = poldegree(P), Q, a, b);
  if(d < 0 || d % 2, return(0)); if(d == 0, return(1));
  Q = P / pollead(P);
  if(d == 2, a = polcoef(Q, 1) / 2; return(Q == ('x + a)^2));
  a = polcoef(Q, 3) / 2; b = (polcoef(Q, 2) - a^2) / 2;
  Q == ('x^2 + a * 'x + b)^2 };
cls_cover(f, c) = {
  my(N = simplify(c[2][1] * 'y^2), R = c[1], F = subst(f, 'x, 'th), P, fa, d);
  if(variable(N) != 'x || poldegree(N, 'y) > 0, error("cls_cover: y in numerator"));
  \\ on each factor K_j of L: 4 N(t) - th R(t) = c_j * square in K_j[t]
  fa = factor(F)[, 1];
  foreach(fa, g, my(Pj = if(poldegree(g) == 1, 4 * N + polcoef(g, 0) / polcoef(g, 1) * R, 4 * N * Mod(1, g) - Mod('th, g) * R));
    if(!iscsquare(Pj), error("cls_cover: not c * square on a factor")));
  P = 4 * N * Mod(1, F) - Mod('th, F) * R;
  forstep(t0 = 0, 200, 1, d = subst(P, 'x, t0); if(Lnorm(d) != 0, break));
  if(Lnorm(d) == 0, error("cls_cover: no invertible value"));
  if(!issquare(Lnorm(d)), error("cls_cover: norm not a square"));
  d };
\\ delta of a point P' = [X, Y] on y^2 = f(x)
cls_point(f, P) = {
  my(F = Mod(1, subst(f, 'x, 'th)), X);
  if(P == [0], return(F));
  X = P[1];
  if(subst(f, 'x, X) != 0, return(X - Mod('th, F.mod)));
  \\ 2-torsion: on the factor (th - X) use g(X), g = f / (x - X); elsewhere X - th
  my(g = f \ ('x - X), fa = factor(subst(f, 'x, 'th))[, 1], mods = vector(#fa, j, if(subst(fa[j], 'th, X) == 0, Mod(subst(g, 'x, X), fa[j]), Mod(X - 'th, fa[j]))));
  my(z = mods[1]); for(j = 2, #mods, z = chinese(z, mods[j])); Mod(lift(z), F.mod) };

\\ characters: n pairs [l, r], l prime >= l0, f(r) = 0 mod l, l not dividing disc(f) * extra
mkchars(f, n, l0, extra) = {
  my(out = List(), D = poldisc(f) * extra, l = nextprime(l0));
  while(#out < n, if(D % l, foreach(polrootsmod(f, l), r, if(#out < n, listput(out, [l, lift(r)])))); l = nextprime(l + 1));
  Vec(out) };
\\ F2 vector of delta (0 = square at that character); returns 0 if delta vanishes at a character
chvec(CH, d) = {
  my(a = lift(d), v = vector(#CH));
  for(i = 1, #CH, my(l = CH[i][1], z = subst(a, 'th, Mod(CH[i][2], l)));
    if(z == 0, return(0)); v[i] = (1 - kronecker(lift(z), l)) / 2);
  v };

\\ ---------- the algebra L with its factors, and reduction of classes modulo squares ----------
\\ A = [f, F (= f in th), comps, S]: comps[j] = [fa_j (in th), nf_j or 0 (linear), [nf, Mod(a, P)] map, back map, prime ideals]
\\ S: rational primes allowed in reductions (2, 3, bad primes, primes of the basis norms)
mkalg(f, S) = {
  my(F = subst(f, 'x, 'th), fa = factor(F)[, 1], comps = vector(#fa), D = abs(poldisc(F)));
  foreach(S, p, D /= p^valuation(D, p));
  if(D != 1, error("mkalg: disc(f) has primes outside S: ", D));
  for(j = 1, #fa,
    if(poldegree(fa[j]) == 1, comps[j] = [fa[j], 0, 0, 0, []]; next);
    my(R = nfinit([fa[j], S], 3), nf = R[1], back = modreverse(R[2]));
    comps[j] = [fa[j], nf, R[2], back, []]);
  [f, F, comps, S] };
\\ register the prime ideals above the primes of S in every component (call after S is final)
alg_primes(A) = {
  for(j = 1, #A[3], my(c = A[3][j]); if(c[2] == 0, next);
    A[3][j][5] = concat(apply(p -> idealprimedec(c[2], p), A[4])));
  A };
\\ component j of an L-element, as an element of nf_j (polmod in nf.pol) or a rational
comp(A, d, j) = { my(c = A[3][j], a = lift(d) % c[1]);
  if(c[2] == 0, return(polcoef(a, 0, 'th)));   \\ fa monic linear: a is constant
  subst(a, 'th, c[3]) };
\\ back to L by CRT
fromcomps(A, z) = {
  my(m = vector(#z, j, my(c = A[3][j]);
    if(c[2] == 0, Mod(z[j], c[1]), Mod(subst(lift(z[j]), variable(c[2].pol), c[4]), c[1]))));
  my(r = m[1]); for(j = 2, #m, r = chinese(r, m[j]));
  Mod(lift(r), A[2]) };
\\ reduce one component: rational -> squarefree core; field -> delta * alpha^2, alpha in b^{-1} minimising
\\ sum_sigma |sigma(delta)| |sigma(alpha)|^2 (balances the unit part too)
redcomp(A, j, z) = {
  my(c = A[3][j]);
  if(c[2] == 0, return(core(numerator(z) * denominator(z))));
  my(nf = c[2], PR = c[5], n = poldegree(nf.pol), e, bi, Bm, emb, w, G, U, al, r1 = nf.sign[1]);
  e = vector(#PR, k, idealval(nf, z, PR[k]));
  my(Nrest = abs(nfeltnorm(nf, z)) / prod(k = 1, #PR, idealnorm(nf, PR[k])^e[k]));
  if(Nrest != 1, error("redcomp: element has primes outside S, norm part ", Nrest));
  bi = idealfactorback(nf, PR, vector(#PR, k, -(e[k] \ 2)));   \\ b^{-1}, b = prod pr^floor(e/2)
  Bm = nf.zk * bi;                                            \\ Z-basis of b^{-1} as polynomials
  localprec(max(100, 3 * ceil(abs(log(abs(nfeltnorm(nf, z)) + 2)) / log(10)) + 60));
  w = abs(nfeltembed(nf, z));
  emb = vector(n, i, nfeltembed(nf, Bm[i]));
  G = matrix(n, n, i, k, sum(s = 1, #w, (if(s <= r1, 1, 2)) * w[s] * real(emb[i][s] * conj(emb[k][s]))));
  U = qflllgram(G);
  my(best = 0, bn = 0);
  for(i = 1, n, my(a = Bm * U[, i], zz = nfeltmul(nf, z, nfeltpow(nf, a, 2)), sz = vecsum(abs(nfeltembed(nf, zz))));
    if(best == 0 || sz < bn, best = zz; bn = sz));
  if(type(best) == "t_COL", best = nfbasistoalg(nf, best));
  best };
redcls(A, d) = fromcomps(A, vector(#A[3], j, redcomp(A, j, comp(A, d, j))));

\\ ---------- cover of a class ----------
\\ matrices of Q_k(u) = coefficient of th^k in delta (u0 + u1 th + u2 th^2)^2
qmats(A, d) = {
  my(M = vector(3, k, matrix(3, 3)));
  for(a = 0, 2, for(b = 0, 2, my(z = lift(d * Mod('th^(a + b), A[2])));
    for(k = 0, 2, M[k + 1][a + 1, b + 1] = polcoef(z, k, 'th))));
  M };
\\ factor |n| using the primes of S first; the cofactor is factored by PARI (should be small)
sfactor(n, S) = {
  n = abs(n); my(P = List(), E = List());
  foreach(S, p, my(v = valuation(n, p)); if(v, listput(P, p); listput(E, v); n /= p^v));
  if(n > 1, my(F = factor(n)); for(i = 1, #F~, listput(P, F[i, 1]); listput(E, F[i, 2])));
  if(#P == 0, return(matrix(0, 2)));
  matconcat([Col(Vec(P)), Col(Vec(E))]) };
\\ quartic cover of class d (reduced) on the member with minimal model E (c4 = E.c4).
\\ Returns [G, Mt, Mp, M0, M1] (old (s,t)~ = Mt (m,n)~; u = Mp [s^2, s t, t^2]~), or 0 if the conic has no point.
cover3(A, d, c4) = {
  my(M = qmats(A, d), M2 = M[3], den = denominator(content(M2)), G2, sol, Mp, U, Q, c, b, r);
  G2 = M2 * den; G2 /= content(G2);
  sol = qfsolve([G2, sfactor(matdet(G2), A[4])]);
  if(type(sol) != "t_COL", return(0));
  Mp = qfparam(G2, sol);
  U = Mp * ['x^2, 'x, 1]~;
  Q = -(U~ * M[2] * U);
  c = content(Q); Q /= c;                       \\ c may be a fraction n/d: c Q = n d Q / d^2
  b = core(numerator(c) * denominator(c), 1); Q *= b[1];
  r = qminred(Q, c4);
  [r[1], r[2], Mp, M[1], M[2]] };
\\ exact check of (m, n) on a kind-3 cover; data = [Mt, Mp, M0, M1, v, f]; returns the point on the member or 0
coverpoint3(cv, m, n) = {
  my(G = cv[3], g, dt = cv[4], t, u, q1, X, Y2, Y);
  g = sum(k = 0, 4, polcoef(G, k) * m^k * n^(4 - k));
  if(!issquare(g), return(0));
  t = dt[1] * [m, n]~;
  u = dt[2] * [t[1]^2, t[1] * t[2], t[2]^2]~;
  q1 = u~ * dt[4] * u;
  if(q1 == 0, return(0));
  X = -(u~ * dt[3] * u) / q1;
  Y2 = subst(dt[6], 'x, X);
  if(!issquare(Y2, &Y), error("coverpoint3: non-point m=", m, " n=", n));
  ellchangepointinv([X, Y], dt[5]) };

\\ ---------- all cosets of Sel2(E)/delta(E(Q)_tors) for one member ----------
\\ E: minimal model of the member; KP: known points ON E (may be []); i: member index (for the cover records).
\\ Returns [covers, info]: covers = list of [3, i, G, [Mt, Mp, M0, M1, v, f], x, new] (x: coset label over the kept
\\ ell2cover basis indices, new = 1 if the class is outside span(torsion, KP)); info = [dimSel2, dimTors, dimKnown,
\\ #cosets, #new, keep (basis indices), seconds].
f2rank(V) = if(#V == 0, 0, matrank(Mod(Mat(V~)~, 2)));
memcovers(E, KP, i, verbose) = {
  my(t0 = getwalltime(), F = fmodel(E), v = F[1], f = F[2], C = ell2cover(E), D, S, A, CH, VD, T, VT, VK, keep = List(), span, nb, out = List(), dS, dT, dK, nnew = 0);
  if(#C == 0, return([[], [0, 0, 0, 0, 0, [], 0]]));
  D = vector(#C, j, cls_cover(f, C[j]));
  S = setunion([2, 3], factor(abs(E.disc))[, 1]~);
  foreach(D, d, S = setunion(S, factor(abs(Lnorm(d)))[, 1]~));
  S = select(p -> p > 1, S);
  A = alg_primes(mkalg(f, S));
  T = elltors(E)[3];
  my(TC = vector(#T, j, cls_point(f, ellchangepoint(T[j], v))), KC = vector(#KP, j, cls_point(f, ellchangepoint(KP[j], v))));
  my(ok = 0, tries = 0);
  while(!ok, tries++; if(tries > 5, error("memcovers: characters vanish"));
    CH = mkchars(f, 96, 10^6 * tries + random(10^5), 1);
    VD = vector(#D, j, chvec(CH, D[j])); VT = vector(#TC, j, chvec(CH, TC[j])); VK = vector(#KC, j, chvec(CH, KC[j]));
    ok = (#select(z -> type(z) != "t_VEC", concat(concat(VD, VT), VK)) == 0));
  dS = f2rank(VD);
  if(dS < #D, error("memcovers: ell2cover classes dependent"));
  if(f2rank(concat(VD, VT)) > dS, error("memcovers: a torsion class is outside the ell2cover span"));
  if(f2rank(concat(VD, VK)) > dS, error("memcovers: a known point's class is outside the ell2cover span"));
  dT = f2rank(VT);
  span = VT;
  for(j = 1, #D, my(ns = concat(span, [VD[j]])); if(f2rank(ns) > f2rank(span), span = ns; listput(keep, j)));
  keep = Vec(keep); nb = #keep;
  my(base = concat(VT, VK)); dK = f2rank(base) - dT;
  forvec(x = vector(nb, j, [0, 1]), if(x == vector(nb), next);
    my(d = Mod(1, A[2]), w = vector(#CH), isnew, r, cv);
    for(k = 1, nb, if(x[k], d *= D[keep[k]]; w = (w + VD[keep[k]]) % 2));
    isnew = f2rank(concat(base, [w])) > f2rank(base);
    r = redcls(A, d);
    cv = cover3(A, r, E.c4);
    if(cv == 0, error("memcovers: conic insoluble for a Selmer class (bug)"));
    if(!issquare(Lnorm(r)), error("memcovers: reduced class norm not a square"));
    nnew += isnew;
    listput(out, [3, i, cv[1], [cv[2], cv[3], cv[4], cv[5], v, f], x, isnew]));
  my(info = [dS, dT, dK, #out, nnew, keep, (getwalltime() - t0) / 1000.]);
  if(verbose, print("# member ", i, " tors ", elltors(E)[2], ": Sel2 dim ", dS, " (ell2cover basis), torsion dim ", dT,
                    ", known points add ", dK, ": ", #out, " cosets mod torsion, ", nnew, " outside span(torsion, known)  [", Strprintf("%.1f", info[7]), " s]"));
  [Vec(out), info] };

\\ point on the input curve M -> member i (minimal model); CL from members() (needs the 4th entry, pts.gp >= this commit)
tomember(CL, i, P) = { my(Q = ellchangepoint(P, CL[1])); if(i > 1, Q = ellisogenyapply(CL[2][i][4], Q)); ellchangepoint(Q, CL[2][i][2]) };
\\ all members: [CV, infos]; KPM = known points on M
covers_all(CL, KPM, verbose) = {
  my(out = [], infos = vector(#CL[2]));
  for(i = 1, #CL[2], my(E = CL[2][i][1], KP = vector(#KPM, j, tomember(CL, i, KPM[j])), R = memcovers(E, KP, i, verbose));
    out = concat(out, R[1]); infos[i] = R[2]);
  [out, infos] };
