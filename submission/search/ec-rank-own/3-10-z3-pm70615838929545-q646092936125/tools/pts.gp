\\ pts.gp -- host side of the GPU point search on 2-coverings (hunt-points brief, helpers/hunt-points.md).
\\ Own code. Inputs used for testing are calibration-only frozen records; nothing here is ever submitted.
\\
\\ Pipeline for a curve M (minimal model):
\\   CL = members(M)              isogeny class: minimal model of every member + transport data back to M
\\   CV = covers(CL)              2-covering quartics y^2 = G(x) for every member:
\\                                 full 2-torsion member: Sel2 classes labelled, torsion image removed, one quartic per
\\                                 nonzero coset of Sel2/delta(E(Q)_tors) (basis cosets: ell2cover quartic; sums:
\\                                 conic parametrisation + minimal model + reduction);
\\                                 otherwise: the ell2cover basis quartics only.
\\   writequartics(CV, file)      "qid a4 a3 a2 a1 a0" (G(m,n) = sum a_k m^k n^(4-k)) for the GPU sieve
\\   coverpoint(CV[q], m, n)      exact check of a hit + map to the member (0 if not a point)
\\   topoint(CL, CV[q], P)        transport to M
\\ Run gp from this directory with -D parisizemax=2000000000.

if(type(PRUNE) == "t_POL", PRUNE = 0);
read("search/selmer/oracle.gp");
read("search/selmer/cassels.gp");
read("search/points/allcov.gp");
toshort(E) = [1, -E.b2/12, -E.a1/2, -E.a3/2 + E.a1*E.b2/24];

\\ ---------- isogeny class ----------
\\ [v, [[Emi, w, dual], ...]]: Emi minimal model of member i, w the change Ei(short, from ellisomat) -> Emi,
\\ dual the isogeny Ei -> Ms (0 for i = 1, Ms = M in short form).
members(M) = {
  my(v = toshort(M), Ms = ellchangecurve(M, v), L = ellisomat(Ms)[1], out = vector(#L));
  for(i = 1, #L, my(Ei = ellinit(L[i][1]), w, Emi = ellminimalmodel(Ei, &w));
    out[i] = [Emi, w, if(i == 1, 0, L[i][3]), L[i][2]]);   \\ [4]: isogeny Ms -> Ei (allcov pruning)
  [v, out, M] };
transport(CL, i, P) = {
  my(mb = CL[2][i], Q = ellchangepointinv(P, mb[2]));
  if(i > 1, Q = ellisogenyapply(mb[3], Q));
  ellchangepointinv(Q, CL[1]) };

\\ ---------- full 2-torsion model y^2 = (x-e1)(x-e2)(x-e3) ----------
emodel(E) = {
  my(v = [1/2, 0, -E.a1/2, -E.a3/2], E2 = ellchangecurve(E, v),
     r = nfroots(, 'x^3 + E.b2*'x^2 + 8*E.b4*'x + 16*E.b6));
  if(#r < 3, return(0));
  if(E2[1..5] != [0, E.b2, 0, 8*E.b4, 16*E.b6], error("emodel"));
  [E2, v, vecsort(r)] };
sqf(a) = core(numerator(a) * denominator(a));
\\ Kummer map on y^2 = prod(x - e_i)
delta(e, P) = {
  if(P == [0], return([1, 1, 1]));
  vector(3, i, my(a = P[1] - e[i]);
    if(a == 0, my(j = (i % 3) + 1, k = ((i + 1) % 3) + 1); a = (e[i] - e[j]) * (e[i] - e[k]));
    sqf(a)) };
clsmul(a, b) = vector(3, i, core(a[i] * b[i]));

\\ class of an ell2cover cover [R, P] on E (min model) in the e-model (E2, v, e): x_E = N(x)/R(x)
covclass(EM, c) = {
  my(v = EM[2], u = v[1], r = v[2], e = EM[3], N = simplify(c[2][1] * 'y^2), R = c[1], d = vector(3));
  if(variable(N) != 'x || poldegree(N, 'y) > 0, error("cover map has y in numerator"));
  for(i = 1, 3, my(Ni = N - (r + u^2 * e[i]) * R, lc = pollead(Ni));
    if(!issquare(Ni / lc), error("covclass: not c*square, i=", i));
    d[i] = sqf(lc));
  d };

\\ F2 coordinates of classes over a prime list S (with -1); class = (d1, d2, d3), d3 = d1 d2 mod squares
f2vec(S, d) = { concat(vector(#S, j, if(S[j] == -1, d[1] < 0, valuation(d[1], S[j]) % 2)),
                       vector(#S, j, if(S[j] == -1, d[2] < 0, valuation(d[2], S[j]) % 2))) };

\\ ---------- quartic from a class via the conic ----------
\\ x - e1 = d1 u1^2/u0^2, x - e2 = d2 u2^2/u0^2  (conic C: d1 u1^2 - d2 u2^2 = (e2 - e1) u0^2),
\\ x - e3 = d3 * square  <=>  W^2 = d3 (d1 U1^2 + (e1 - e3) U0^2)  on the parametrisation (U0,U1,U2) = Mp [t^2,t,1].
\\ Returns [G (minimal, reduced quartic in x), Mt (2x2: (t1,t2)~ = Mt (m,n)~), Mp, e, d] or 0 when the conic has no point.
\\ c4 = c4 of the member's minimal model: the minimal quartic has I(G) = c4.
conicquartic(e, d, c4) = {
  my(G = matdiagonal([-(e[2] - e[1]), d[1], -d[2]]), F = factor(abs((e[2] - e[1]) * d[1] * d[2])), sol, Mp, U, Q, c, b, r);
  sol = qfsolve([G, F]);
  if(type(sol) != "t_COL", return(0));
  Mp = qfparam(G, sol);
  U = Mp * ['x^2, 'x, 1]~;
  Q = d[3] * (d[1] * U[2]^2 + (e[1] - e[3]) * U[1]^2);
  c = content(Q); Q /= c;                     \\ Q primitive now, true quartic = c * Q
  b = core(c, 1); Q *= b[1];                  \\ c = b[1] * b[2]^2: W' = W / b[2]
  r = qminred(Q, c4);
  [r[1], r[2], Mp, e, d] };

\\ ---------- minimisation and reduction of y^2 = G(x) ----------
qI(G) = { my(a=polcoef(G,4), b=polcoef(G,3), c=polcoef(G,2), d=polcoef(G,1), e=polcoef(G,0)); 12*a*e - 3*b*d + c^2 };
qJ(G) = { my(a=polcoef(G,4), b=polcoef(G,3), c=polcoef(G,2), d=polcoef(G,1), e=polcoef(G,0)); 72*a*c*e + 9*b*c*d - 27*a*d^2 - 27*e*b^2 - 2*c^3 };
\\ G(a x + b, c x + d) (homogeneous substitution, dehomogenised at y = 1)
qsub(G, A) = { sum(k = 0, 4, polcoef(G, k) * (A[1,1]*'x + A[1,2])^k * (A[2,1]*'x + A[2,2])^(4 - k)) };
\\ one p-adic step: index-D sublattice A (Hermite form, D = p^j) with G o A == 0 mod p^2 D^2; I drops by p^4
qminstep(G, p) = {
  for(j = 1, 3, my(D = p^j);
    fordiv(D, a, my(dd = D / a); for(b = 0, a - 1, my(A = [a, b; 0, dd], H = qsub(G, A));
      if(content(H) % (p^2 * D^2) == 0, return([H / (p^2 * D^2), A])))));
  0 };
\\ minimise (I(G) -> c4) and SL2(Z)-reduce; returns [G', Mt] with old x = Mt . x'
qminred(G, c4) = {
  my(Mt = matid(2), ratio, M2, G2);
  ratio = qI(G) / c4;
  if(type(ratio) != "t_INT", error("qminred: I/c4 not integral: ", ratio));
  foreach(factor(ratio)[, 1]~, p,
    if(p > 1000,                                    \\ large p: PARI's minimisation (the residue loop below is O(p))
      my(m, H = hyperellminimalmodel(G, &m, [p]));
      if(type(H) == "t_VEC", if(H[2] != 0, error("hyperellminimalmodel returned a Q part")); H = H[1]);
      if(type(m) == "t_VEC" && #m >= 3 && m[3] != 0, error("hyperellminimalmodel: y-shift"));
      G = H; Mt = Mt * mobmat(m); next);
    while(qI(G) / c4 % p^4 == 0, my(st = qminstep(G, p)); if(st == 0, break); G = st[1]; Mt = Mt * st[2]));
  G2 = hyperellred(G, &M2);
  if(type(G2) == "t_VEC", if(G2[2] != 0, error("hyperellred returned a Q part")); G2 = G2[1]);
  Mt = Mt * mobmat(M2);
  [G2, Mt] };
\\ hyperellchangecurve convention: old x = (a x' + b)/(c x' + d); M = 0 / scalar means identity
mobmat(m) = if(type(m) == "t_VEC" && #m >= 2 && type(m[2]) == "t_MAT", m[2], matid(2));
\\ Stoll-Cremona covariant point z(G) in the upper half plane (weights t_j minimising a^2 disc(Q_t)^2 / prod t_j)
sc_z(G) = {
  localprec(80);
  my(r = polroots(G), n = #r, t = vector(n, j, 1.), A, B, C, D, g, tn);
  for(it = 1, 500,
    A = vecsum(t); B = -2 * sum(j = 1, n, t[j] * real(r[j])); C = sum(j = 1, n, t[j] * abs(r[j])^2); D = 4*A*C - B^2;
    g = vector(n, j, 4 * (C + A * abs(r[j])^2 + B * real(r[j])));
    tn = vector(n, j, D / ((n / 2) * g[j])); tn /= vecsum(tn);
    tn = vector(n, j, sqrt(tn[j] * t[j])); tn /= vecsum(tn);
    if(normlp(tn - t) < 1e-30, t = tn; break); t = tn);
  A = vecsum(t); B = -2 * sum(j = 1, n, t[j] * real(r[j])); C = sum(j = 1, n, t[j] * abs(r[j])^2);
  (-B + I * sqrt(4*A*C - B^2)) / (2*A) };
\\ search box aspect: points with small covariant |m - z n|^2 / Im z lie in |m| <= B s, 1 <= n <= B / s, s = sqrt(Im z)
qaspect(G) = { if(poldegree(G) < 4, return([1., 0.])); my(z = sc_z(G)); [sqrt(imag(z)), real(z)] };

\\ ---------- covers of a class ----------
\\ cover record: [kind, member, G, data]
\\   kind 1: data = [R, P] (ell2cover)        kind 2: data = [Mt, Mp, e, d, EM]
\\ Cassels-Tate kernel classes for a full-2-torsion member (hunt-selmer's cassels.gp, validated against ellrank):
\\ returns the F2 kernel of the pairing as a list of triples spanning it, or 0 on failure.
casskernel(e) = {
  my(r = e[2] - e[1], s = e[3] - e[1], R, C, BAS, K, out = List());
  for(attempt = 1, 4,
    R = iferr(casselsmat(r, s), err, 0);
    if(R != 0 && R[1] == R[1]~ && sum(i = 1, #R[1], R[1][i, i]) == 0, break); R = 0);
  if(R == 0, return(0));
  C = R[1]; BAS = R[3];
  K = lift(matker(Mod(C, 2)));
  for(c = 1, #K, my(d = [1, 1, 1]); for(j = 1, #BAS, if(K[j, c], d = clsmul(d, BAS[j]))); listput(out, d));
  Vec(out) };
\\ ALLCOV = 1: members without full rational 2-torsion get one quartic per coset of Sel2/delta(tors) (allcov.gp),
\\ classes outside span(torsion, KNOWN) flagged (KNOWN = points on the input curve, default []).
if(type(ALLCOV) == "t_POL", ALLCOV = 0);
if(type(KNOWN) == "t_POL", KNOWN = []);
allcov_member(CL, i, verbose) = { my(KP = vector(#KNOWN, j, tomember(CL, i, KNOWN[j]))); memcovers(CL[2][i][1], KP, i, verbose)[1] };
covers(CL, verbose, ml) = {
  if(type(ml) != "t_VEC", ml = [1..#CL[2]]);    \\ ml: members to build (default all)
  if(PRUNE, return(covers_pruned(CL, verbose, ml)));
  my(out = List());
  foreach(ml, i, my(E = CL[2][i][1], EM = emodel(E), C);
    if(EM == 0 && ALLCOV, foreach(allcov_member(CL, i, verbose), c, listput(out, c)); next);
    C = ell2cover(E);
    if(EM == 0 || #C == 0,
      foreach(C, c, listput(out, [1, i, c[1], c]));
      if(verbose, print("# member ", i, " tors ", elltors(E)[2], ": not full 2-torsion, ", #C, " basis covers"));
      next);
    my(e = EM[3], cl = vector(#C, j, covclass(EM, C[j])), T = elltors(E)[3], tc, S, V, TV, keep = List(), span, nb);
    tc = vector(#T, j, delta(e, ellchangepoint(T[j], EM[2])));
    S = Set(concat([-1], concat(vector(#cl, j, concat(vector(3, k, factor(abs(cl[j][k]))[, 1]~))))));
    S = Set(concat(S, concat(vector(#tc, j, concat(vector(3, k, factor(abs(tc[j][k]))[, 1]~))))));
    S = Set(concat([-1], select(p -> p > 1, S)));
    TV = vector(#tc, j, f2vec(S, tc[j]));
    span = TV;
    for(j = 1, #cl, my(w = f2vec(S, cl[j]), ns = concat(span, [w]));
      if(matrank(Mod(Mat(ns~)~, 2)) > matrank(Mod(Mat(span~)~, 2)) || #span == 0, span = ns; listput(keep, j)));
    nb = #keep;
    if(verbose, print("# member ", i, " tors ", elltors(E)[2], ": Sel2 basis ", #C, ", torsion rank ", matrank(Mod(Mat(TV~)~, 2)), ", cosets 2^", nb, "-1"));
    forvec(x = vector(nb, j, [0, 1]), if(x == vector(nb), next);
      my(js = select(k -> x[k], [1..nb]));
      if(#js == 1,
        my(c = C[keep[js[1]]]); listput(out, [1, i, c[1], c]),
        my(d = [1, 1, 1], q);
        foreach(js, k, d = clsmul(d, cl[keep[k]]));
        q = conicquartic(e, d, E.c4);
        if(q == 0, if(verbose, print("# member ", i, " class ", d, ": conic insoluble (not in Sel2?)")); next);
        listput(out, [2, i, q[1], [q[2], q[3], q[4], q[5], EM]]))));
  Vec(out) };

\\ pruned: full-2-torsion members get one conic quartic per nonzero coset of (Cassels kernel)/delta(torsion);
\\ other members get the ell2cover basis (no pruning possible there yet).
covers_pruned(CL, verbose, ml) = {
  if(type(ml) != "t_VEC", ml = [1..#CL[2]]);
  my(out = List());
  foreach(ml, i, my(E = CL[2][i][1], EM = emodel(E));
    if(EM == 0 && ALLCOV, foreach(allcov_member(CL, i, verbose), c, listput(out, c)); next);
    if(EM == 0,
      my(C = ell2cover(E)); foreach(C, c, listput(out, [1, i, c[1], c]));
      if(verbose, print("# member ", i, " tors ", elltors(E)[2], ": not full 2-torsion, ", #C, " basis covers"));
      next);
    my(e = EM[3], K = casskernel(e), T = elltors(E)[3], tc, S, gens, span, keep = List(), nb);
    if(K == 0, error("casskernel failed on member ", i));
    tc = vector(#T, j, delta(e, ellchangepoint(T[j], EM[2])));
    gens = concat(tc, K);
    S = Set(concat([-1], concat(vector(#gens, j, concat(vector(3, k, factor(abs(gens[j][k]))[, 1]~))))));
    S = Set(concat([-1], select(p -> p > 1, S)));
    span = vector(#tc, j, f2vec(S, tc[j]));
    for(j = 1, #K, my(ns = concat(span, [f2vec(S, K[j])]));
      if(matrank(Mod(Mat(ns~)~, 2)) > matrank(Mod(Mat(span~)~, 2)), span = ns; listput(keep, K[j])));
    nb = #keep;
    if(verbose, print("# member ", i, " tors ", elltors(E)[2], ": Cassels kernel dim ", #K, ", cosets mod torsion 2^", nb, "-1"));
    forvec(x = vector(nb, j, [0, 1]), if(x == vector(nb), next);
      my(d = [1, 1, 1], q);
      for(k = 1, nb, if(x[k], d = clsmul(d, keep[k])));
      q = conicquartic(e, d, E.c4);
      if(q == 0, error("conic insoluble for a Cassels-kernel class ", d));
      listput(out, [2, i, q[1], [q[2], q[3], q[4], q[5], EM]])));
  Vec(out) };

\\ ---------- exact check of a hit (m, n) on cover cv; returns the point on the member's minimal model or 0 ----------
coverpoint(cv, m, n) = {
  if(cv[1] == 3, return(coverpoint3(cv, m, n)));   \\ allcov.gp
  my(G = cv[3], g, w);
  g = sum(k = 0, 4, polcoef(G, k) * m^k * n^(4 - k));
  if(!issquare(g, &w), return(0));
  if(cv[1] == 1,
    if(n == 0, return(0));
    my(P = subst(subst(cv[4][2], 'x, m / n), 'y, w / n^2));
    return(P));
  my(Mt = cv[4][1], Mp = cv[4][2], e = cv[4][3], d = cv[4][4], EM = cv[4][5], t = Mt * [m, n]~, u, X, Y2, Y);
  u = Mp * [t[1]^2, t[1] * t[2], t[2]^2]~;
  if(u[1] == 0, return(0));
  X = e[1] + d[1] * u[2]^2 / u[1]^2;
  Y2 = (X - e[1]) * (X - e[2]) * (X - e[3]);
  if(!issquare(Y2, &Y), error("coverpoint: conic map gave a non-point m=", m, " n=", n));
  ellchangepointinv([X, Y], EM[2]) };

topoint(CL, cv, P) = transport(CL, cv[2], P);

\\ real intervals of x where G(x) >= 0, as strings for the sieve ("inf" for unbounded)
qintervals(G) = {
  my(r = polrootsreal(G), lc = sign(pollead(G)), dg = poldegree(G), out = List(), s, a);
  s = if(dg % 2, -lc, lc);            \\ sign at -infinity
  a = "-inf";
  for(i = 1, #r, my(x = Strprintf("%.17g", r[i]));
    if(s > 0, listput(out, [a, x]));
    s = -s; a = x);
  if(s > 0, listput(out, [a, "inf"]));
  Vec(out) };
\\ one line per cover: qid mmax nmax a4 a3 a2 a1 a0 nint lo hi ...   box: |m| <= B s, 1 <= n <= B / s (square if SQ)
writequartics(CV, file, B, SQ, qs) = {
  if(type(qs) != "t_VEC", qs = [1..#CV]);
  foreach(qs, q, my(G = CV[q][3], as = if(SQ, [1., 0.], qaspect(G)), s = as[1], mm, nn, iv = qintervals(G), str);
    mm = ceil(B * s); nn = max(1, floor(B / s));
    if(!SQ && polcoef(G, 3) == 0 && polcoef(G, 1) == 0, mm = -mm);   \\ even quartic: x -> -x is an automorphism, m >= 0 suffices
    str = Str(q, " ", mm, " ", nn, " ", polcoef(G, 4), " ", polcoef(G, 3), " ", polcoef(G, 2), " ", polcoef(G, 1), " ", polcoef(G, 0), " ", #iv);
    foreach(iv, v, str = Str(str, " ", v[1], " ", v[2]));
    write(file, str)) };
\\ read sieve hits "qid m n" and return the exact points as [qid, m, n, P (member), Q (on M)]; adds n = 0 points
readhits(CL, CV, file, qs) = {
  if(type(qs) != "t_VEC", qs = [1..#CV]);
  my(L = if(file == "", [], readstr(file)), out = List(), bad = 0);
  foreach(L, s, my(v = strsplit(s, " "), q = eval(v[1]), m = eval(v[2]), n = eval(v[3]), P = coverpoint(CV[q], m, n));
    if(P == 0, bad++; next);
    if(!ellisoncurve(CL[2][CV[q][2]][1], P), error("readhits: off member curve q=", q));
    my(Q = topoint(CL, CV[q], P));
    if(!ellisoncurve(CL[3], Q), error("readhits: transported point off the input curve q=", q));
    listput(out, [q, m, n, P, Q]));
  foreach(qs, q, my(P = coverpoint(CV[q], 1, 0)); if(P != 0, listput(out, [q, 1, 0, P, topoint(CL, CV[q], P)])));
  [Vec(out), bad] };

\\ independent subset (greedy on the height Gram determinant; relative threshold)
indep(E, V) = {
  my(S = [], d0 = 1);
  foreach(V, P, if(P == [0] || ellorder(E, P) > 0, next);
    my(T = concat(S, [P]), dd = matdet(ellheightmatrix(E, T)));
    if(dd > 1e-20 * d0 * ellheight(E, P), S = T; d0 = dd));
  S };
