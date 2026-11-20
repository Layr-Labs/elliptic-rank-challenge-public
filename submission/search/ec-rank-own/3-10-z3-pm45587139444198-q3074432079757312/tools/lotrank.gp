\\ lotrank.gp -- build-lottery ranking stage (own code).  Driver sets INFILE (survivor lines "fam p q s300 s2000 s1e4 H=..."),
\\ OUTFILE, T1 (s, ellrank effort 1), T3 (s, effort 3), TARGET (11).  Per curve: minimal model, torsion check, |a_p| fingerprint
\\ against the three public rank-10 test curves (isogeny- and twist-invariant: any match is dropped, never ranked further),
\\ root number, S(1e5), ellrank(E,1) with alarm T1; effort 3 (seeded with the effort-1 points) when the upper bound reaches
\\ TARGET (or 10 when the class minimum H is below the board) and the lower bound does not; for lower >= 10 the isogeny-class
\\ member of smallest H with the same torsion.  Output line: fam p q s300 s2000 s1e4 H w r_lo r_up eff secs S1e5 Hmin_class
\\ tag ainv [..] pts [[x,y],..]
read("search/lottery/fam3.gp");
TORS = Map(["z6k", 6; "z5", 5; "leg", 8; "z3", 3; "z4", 4]);
\\ board heights to beat per (family, rank) (logs/board_parent_verify.json: min(baseline, board_best)); rank >= 11 = new cell
BOARD = Map(["z6k/8", 15.260014; "z6k/9", 18.445122; "z6k/10", 21.987758; "z5/8", 15.726734; "z5/9", 20.436404; "z5/10", 26.003803; "leg/8", 15.323158; "leg/9", 19.016457; "leg/10", 24.92637; "z3/8", 8.064684; "z3/9", 9.512795; "z3/10", 12.159425; "z3/11", 13.999595; "z3/12", 17.864161; "z3/13", 20.199607; "z3/14", 22.493691; "z3/15", 28.255216; "z4/7", 8.339871; "z4/8", 10.189149; "z4/9", 12.166195; "z4/10", 15.3894; "z4/11", 18.746574; "z4/12", 23.595951; "z4/13", 32.024466; "z6k/7", 12.699047; "z5/7", 13.199808; "leg/7", 12.948901]);
\\ first rank with no frozen cell per family (a new cell scores at any height)
NEWC = Map(["z6k", 11; "z5", 11; "leg", 11; "z3", 16; "z4", 14]);
cellok(f, k, h) = if(k >= mapget(NEWC, f), 1, if(k < 7 || !mapisdefined(BOARD, Str(f, "/", k)), 0, h < mapget(BOARD, Str(f, "/", k))));
reach(f, lo, up, h) = { for(k = max(lo + 1, 7), up, if(cellok(f, k, h), return(1))); 0 };
afp(E) = { my(a = ellan(E, 1000)); apply(p -> abs(a[p]), primes([5, 1000])) };
TESTFP = apply(v -> afp(ellminimalmodel(famE(v[1], v[2], v[3]))), [["z6k", 686288691, 12859057120], ["z5", -121186667075, 43349646368], ["leg", 33661279586, 5418397845]]);
classminH(M, n) = { my(L = ellisomat(M), h = heightE(M));
  if(type(L) == "t_VEC", for(j = 1, #L[1], my(E2 = ellminimalmodel(ellinit(L[1][j][1]))); if(elltors(E2)[1] == n, h = min(h, heightE(E2))))); h };
fmtpts(P) = Str("[", strjoin(apply(x -> Str("[", x[1], ",", x[2], "]"), P), ","), "]");
{
my(L = readstr(INFILE));
foreach(L, s, if(#s == 0, next);
  my(v = select(x -> #x, strsplit(s, " ")), f = v[1], P = eval(v[2]), Q = eval(v[3]), t0 = getwalltime(), M, n, fp, w, s5, r, eff = 1, hmin = 0., tag = "ok");
  M = ellminimalmodel(famE(f, P, Q)); n = elltors(M)[1];
  if(n != mapget(TORS, f), write(OUTFILE, f, " ", P, " ", Q, " ", v[4], " ", v[5], " ", v[6], " ", Strprintf("%.6f", heightE(M)), " tors ", n, " TORSFAIL"); next);
  fp = afp(M);
  if(fp == TESTFP[1] || fp == TESTFP[2] || fp == TESTFP[3], write(OUTFILE, f, " ", P, " ", Q, " TESTCURVE_CLASS dropped (public rank-10 test curve; never ranked)"); next);
  my(reg = v[#v]); if(#reg < 3 || Vecsmall(reg)[1] != 82 || Vecsmall(reg)[2] != 58, reg = "R:?");
  w = ellrootno(M); s5 = Ssums(M, [10^5])[1];
  my(en = "");
  r = iferr(alarm(T1, ellrank(M, 1)), e, en = errname(e); [-1, -1, 0, []]);
  my(hf = heightE(M));
  if(r[1] >= 7 || r[2] >= 7, hmin = classminH(M, n));
  \\ effort 3 when the upper bound reaches a cell this fibre would score in (new cell >= 11, or H below the board of rank k)
  if(r[1] >= 0 && r[2] > r[1] && reach(f, r[1], r[2], hf),
    my(r3 = iferr(alarm(T3, ellrank(M, 3, r[4])), e, en = errname(e); 0)); if(r3 != 0, r = r3; eff = 3, tag = Str("eff3err:", en)));
  if(r[1] < 0, tag = Str("err:", en));   \\ e_ALARM = the T1 alarm, e_STACK = stack overflow (owner 20:55Z)
  if(r[1] >= 7 && hmin == 0., hmin = classminH(M, n));
  write(OUTFILE, f, " ", P, " ", Q, " ", v[4], " ", v[5], " ", v[6], " ", Strprintf("%.6f", heightE(M)), " ", w, " ", r[1], " ", r[2], " ", eff, " ",
        Strprintf("%.2f", (getwalltime()-t0)/1000.), " ", Strprintf("%.2f", s5), " ", Strprintf("%.6f", hmin), " ", tag, " ", reg,
        " ainv [", strjoin(apply(x -> Str(x), Vec(M[1..5])), ","), "] pts ", fmtpts(if(r[1] > 0, r[4], []))));
}
quit
