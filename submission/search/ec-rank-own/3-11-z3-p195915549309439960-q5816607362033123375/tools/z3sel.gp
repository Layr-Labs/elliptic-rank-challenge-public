\\ z3sel.gp -- closed-form 3-isogeny Selmer bound for the Tate Z/3 family  E: y^2 + a x y + b y = x^3,  T = (0,0) of order 3.
\\ Sieve pair (P, Q) = (b, a^3).  z3bound(P, Q) = [bound, s_phi, s_phihat] with bound = s_phi + s_phihat - 1 >= rank E(Q).
\\ Derived from search/selmer/oracle3.gp (same phi, same Kummer map, same Greenberg-Wiles step); see helpers/z3sel.md.
\\
\\ phi: E -> E' = E/<T>; Kummer map of the dual (kernel mu3): R -> y(R) in Q^*/Q^*3 (y = tangent line at T); kappa(T) = 1/b.
\\ Work on the reduced model (no prime l with l | a and l^3 | b; (a,b) ~ (u a, u^3 b)). D = a^3 - 27 b, Delta = b^3 D.
\\ S = {3} u primes(b) u primes(D); Sel^phihat sits in <S> (-1 is a cube). Local image W_l at l != 3:
\\   l | b, l !| a  (split I_{3k})                 : everything          (no condition)
\\   l | D, l !| b  (I_m, T in the identity comp.) : trivial             (d must be an l-adic cube)
\\   l | a and l | b (additive IV / IV*)          : <b>                 (E = E_0 + <T>, E_0 3-divisible)
\\ at l = 3, with al = v3(a), be = v3(b) (reduced: al >= 1 => be <= 2), u3(x) = class of x in Q_3^*/Q_3^*3:
\\   al = 0: be = 0 -> units (good ordinary, flat mu3);  be >= 1 -> everything (Tate curve, q = cube)
\\   al = 1: be = 0 -> <b>, except s = b/(a/3)^3 = 19 mod 27 -> units;  be = 1 -> everything if a/3 = b/3 mod 3, else <b>;
\\           be = 2 -> everything
\\   al >= 2: be = 0 -> <b>;  be = 1 -> everything;  be = 2 -> <b>
\\ s_phihat = #S - rank_F3(conditions);  s_phi = s_phihat + sum_{l in S} (h_l - dim W_l - 1),  h_l = dim_F3 Q_l^*/Q_l^*3.

\\ class of u (l-adic unit, rational) in the unit part of Q_l^*/Q_l^*3:  l = 3: u mod 9 (1,8 -> 0; 2,7 -> 1; 4,5 -> 2),
\\ l = 1 mod 3: cubic residue character relative to the cube root of unity z (fixed per l)
z3_cls9(u) = { my(r = u % 9); [0,1,0,2,2,0,1,0][r] };
z3_chi(u, l, z) = { my(c = Mod(u, l)^((l-1)/3)); if(c == 1, 0, if(c == z, 1, 2)) };

\\ reduce (a,b) -> (a/u, b/u^3) at every prime l with l | a and l^3 | b
z3_reduce(a, b) = {
  my(g = gcd(a, b));
  if(g == 1, return([a, b]));
  foreach(factor(g)[,1]~, l, while(a % l == 0 && b % l^3 == 0, a /= l; b /= l^3));
  [a, b] };

\\ main entry: (P,Q) = (b, a^3) as in the sieve lines
z3bound(P, Q) = { my(a); if(!ispower(Q, 3, &a), error("z3bound: Q is not a cube: ", Q)); z3boundab(a, P) };

\\ local image type at l on the reduced model: 0 trivial, 1 units (l = 3 only), 2 everything, 3 <b> (a line; v_l(b) != 0 mod 3)
z3_type(l, a, b) = {
  if(l != 3,
    if(b % l, return(0));                                   \\ l | D only: I_m with T in the identity component
    return(if(l % 3 == 2 || a % l, 2, 3)));                 \\ split I_{3k}: everything; additive: <b> (everything if l = 2 mod 3)
  my(al = valuation(a, 3), be = valuation(b, 3), W);
  if(al == 0, return(if(be == 0, 1, 2)));
  \\ the 3 | a rows were read off the oracle and validated for v3(a) <= 25 and v3(a^3 - 27b) <= 40: refuse beyond (fail-safe)
  if(al > 25 || valuation(a^3 - 27*b, 3) > 40, error("z3bound: 3-adic depth beyond the validated range: ", [al, valuation(a^3 - 27*b, 3)]));
  if(al == 1, W = if(be == 0, if(lift(Mod(b, 27) / Mod(a/3, 27)^3) == 19, 1, 3), be == 1, if((a/3 - b/3) % 3 == 0, 2, 3), 2),
    W = if(be == 1, 2, 3));
  if(W == 3 && be == 0, W = if(z3_cls9(b) == 0, 0, 1));    \\ <b> for a 3-adic unit b: trivial (b a cube) or units
  W };

\\ per-prime local data [reduced [a,b], S, types] (used by the validation scripts)
z3loc(a0, b0) = {
  my(ab = z3_reduce(a0, b0), a = ab[1], b = ab[2], D = a^3 - 27*b);
  my(S = setunion(setunion(factor(abs(b))[,1]~, factor(abs(D))[,1]~), [3]));
  foreach(S, l, if(l > 2^64 && !isprime(l), error("z3bound: factor not proven prime: ", l)));   \\ BPSW is proven below 2^64
  [ab, S, vector(#S, i, z3_type(S[i], a, b))] };

z3boundab(a0, b0) = {
  if(a0 == 0 || b0 == 0 || a0^3 == 27*b0, error("z3bound: singular or not in the family: ", [a0, b0]));
  my(L = z3loc(a0, b0), a = L[1][1], b = L[1][2], S = L[2], n = #S, rows = List(), corr = 0);
  for(i = 1, n, my(l = S[i], W = L[3][i], h = if(l == 3 || l % 3 == 1, 2, 1), U, vb, cb);
    corr += h - (W == 1 || W == 3) - 2*(W == 2)*(h == 2) - (W == 2)*(h == 1) - 1;
    if(W == 2, next);
    if(h == 1, my(r = vector(n)); r[i] = 1; listput(rows, r); next);     \\ l = 2 mod 3, trivial: v_l(d) = 0 mod 3
    \\ h = 2: classes (v_l, unit class); U = unit classes of the basis primes (0 at l itself)
    if(l == 3, U = vector(n, j, if(j == i, 0, z3_cls9(S[j]))),
      my(z = lift(polrootsmod('x^2 + 'x + 1, l)[1])); U = vector(n, j, if(j == i, 0, z3_chi(S[j], l, z))));
    if(W == 0, my(r = vector(n)); r[i] = 1; listput(rows, r); listput(rows, U); next);
    if(W == 1, my(r = vector(n)); r[i] = 1; listput(rows, r); next);
    \\ W = <b>, class (vb, cb), vb != 0 mod 3: d in <b> iff vb*u(d) - cb*v(d) = 0
    vb = valuation(b, l); cb = if(l == 3, z3_cls9(b / 3^vb), z3_chi(b / l^vb, l, lift(polrootsmod('x^2 + 'x + 1, l)[1])));
    listput(rows, vector(n, j, vb*U[j] - cb*(j == i))));
  my(rk = if(#rows, matrank(Mod(matrix(#rows, n, i, j, rows[i][j]), 3)), 0), sh = n - rk);
  [2*sh + corr - 1, sh + corr, sh] };
