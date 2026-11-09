\\ kub_rank.gp -- rank stage for the candidates printed by `kub sieve` (PARI/GP 2.17.4; own code).
\\ usage: KUB_FAM=6 KUB_IN=cand.txt [KUB_MINRANK=7] [KUB_EFFORT=1] [EC_OUT=found] gp -q kub_rank.gp
\\ input line: "p q S300 S2000 S10000 H"   (t = p/q); output line: "p q S300 S2000 S10000 Hexact lb ub npts seconds"
\\ families: 5 = Z/5, 6 = Z/6, 42 = Z/2xZ/4 (Tate normal forms, see kub.c); the torsion group is checked.
default(parisize, 200000000);
default(parisizemax, 12000000000);
fam = eval(getenv("KUB_FAM"));
fn = getenv("KUB_IN");
minr = if(getenv("KUB_MINRANK"), eval(getenv("KUB_MINRANK")), 99);
effort = if(getenv("KUB_EFFORT"), eval(getenv("KUB_EFFORT")), 1);
outdir = if(getenv("EC_OUT"), getenv("EC_OUT"), "found");
tmax = if(getenv("KUB_TMAX"), eval(getenv("KUB_TMAX")), 300);
\\ stage 4 (cheap screen before ellrank): dS = S(10^5) - S(10^4) = sum over the primes 10^4 < l <= 10^5 of -a_l log(l)/l, 0.07 s per curve.  Measured on Z/5 curves of known rank:
\\ rank 0-2 mean -1.5..0.3 (sd 4-6), rank 6 mean 8.5, rank 7 mean 15.7 (min 6.8), rank 8 15.6, frozen rank-9 records 19.7 (Z/5), 25.1 (Z/6), 24.8 (2xZ/4): a cut dS >= 8 keeps >99.9 % of rank 9
\\ and lets through only ~5 % of the low-rank candidates.  KUB_DSMIN (default off) sets the cut; curves below it are printed with lb = ub = -2 and not ranked.
dsmin = if(getenv("KUB_DSMIN"), eval(getenv("KUB_DSMIN")), -1e9);
\\ optional parity screen: KUB_WSIGN = -1 keeps only curves whose root number is -1 (odd analytic rank, hence odd rank by the parity conjecture) or +1 for even rank; ellrootno costs ~ms.
wsign = if(getenv("KUB_WSIGN"), eval(getenv("KUB_WSIGN")), 0);
famfile = if(getenv("KUB_FAMS"), getenv("KUB_FAMS"), "fams.txt");
LL = readstr(famfile);
FD = 0; AA = vector(5, k, [0]); AW = [1, 2, 3, 4, 6];     \\ blocks for [a1,a2,a3,a4,a6] with weights AW; A4/A6 are optional trailing blocks (family 44 = Legendre Z/2xZ/4)
{
  my(i = 1, s = 0, n = 0);
  while(i <= #LL,
    my(t = strsplit(LL[i], " "));
    if(t[1] == "FAM" && eval(t[2]) == fam,
      FD = eval(t[3]);
      for(k = 1, 3, s = select(x -> #x > 0, strsplit(LL[i + k], " ")); n = eval(s[2]); AA[k] = vector(n, j, eval(s[2 + j])));
      for(k = 4, 5,
        if(i + k <= #LL,
          s = select(x -> #x > 0, strsplit(LL[i + k], " "));
          if(#s > 1 && s[1] == if(k == 4, "A4", "A6"), n = eval(s[2]); AA[k] = vector(n, j, eval(s[2 + j])))));
      break);
    i++);
  if(FD == 0, error("family ", fam, " not found in ", famfile));
}
kmodel(p, q) = vector(5, k, sum(j = 1, #AA[k], AA[k][j] * p^(j-1) * q^(AW[k]*FD - (j-1))));
wantt = if(fam == 4, [4, [4]], if(fam == 5, [5, [5]], if(fam == 6 || fam == 66 || fam == 67, [6, [6]], if(fam == 42 || fam == 44, [8, [4, 2]], if(fam == 8 || fam == 88, [8, [8]], if(fam == 10, [10, [10]], if(fam == 12, [12, [12]], if(fam == 26 || fam == 261 || fam == 265, [12, [6, 2]], [16, [8, 2]]))))))));
cstr(x) = if(denominator(x) == 1, Str(numerator(x)), Str(numerator(x), "/", denominator(x)));
L = readstr(fn);
{
for(ii = 1, #L,
  my(tok, p, q, E, Em, w, tor, t0, r, H, pts, lb, ub, fname, G, dS);
  t0 = getabstime();
  tok = strsplit(L[ii], " "); tok = select(s -> #s > 0, tok);
  if(#tok < 2, next);
  p = eval(tok[1]); q = eval(tok[2]);
  E = ellinit(kmodel(p, q));
  if(E == [], next);
  Em = ellminimalmodel(E, &w);
  tor = elltors(Em); tor = [tor[1], tor[2]];
  if(tor != wantt, print(p, " ", q, " torsion-mismatch ", tor); next);
  H = log(max(abs(Em.c4)^3, Em.c6^2))/12.;
  if(wsign != 0 && ellrootno(Em) != wsign,
    printf("%d %d %s %.5f -3 -3 0 %.2f 0.00\n", p, q, if(#tok > 4, strjoin(tok[3..5], " "), "- - -"), H, (getabstime() - t0)/1000.);
    next);
  dS = 0.;
  if(dsmin > -1e8, forprime(l = 10001, 100000, dS += -ellap(Em, l) * log(l) / l));
  if(dS < dsmin,
    printf("%d %d %s %.5f -2 -2 0 %.2f %.2f\n", p, q, if(#tok > 4, strjoin(tok[3..5], " "), "- - -"), H, (getabstime() - t0)/1000., dS);
    next);
  r = iferr(alarm(tmax, ellrank(Em, effort)), e, [-1, -1, 0, []]);
  if(type(r) == "t_ERROR", r = [-1, -1, 0, []]);        \\ alarm(s, code) RETURNS a t_ERROR when the time limit is hit (it does not throw): without this line r[1] aborted the whole script and every later line of the part file was never ranked
  lb = r[1]; ub = r[2]; pts = r[4];
  printf("%d %d %s %.5f %d %d %d %.2f %.2f\n", p, q, if(#tok > 4, strjoin(tok[3..5], " "), "- - -"), H, lb, ub, #pts, (getabstime() - t0)/1000., dS);
  if(lb >= minr,
    pts = ellsaturation(Em, pts, 100);
    G = ellheightmatrix(Em, pts);
    fname = Str(outdir, "/", fam, "-", #pts, "-", p, "_", q, ".ec");
    write(fname, "E ", strjoin(apply(s -> Str(s), Em[1..5]), " "));
    for(m = 1, #pts, write(fname, "P ", cstr(pts[m][1]), " ", cstr(pts[m][2])));
    write(fname, "#CLAIM rank=", #pts);
    print("FOUND ", fname, " H=", H, " lb=", lb, " npts=", #pts, " det=", matdet(G));
  );
);
}
quit;
