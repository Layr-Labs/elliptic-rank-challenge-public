\\ ptsearch.gp -- one curve: covers of every class member -> GPU sieve -> exact check -> points on the input curve.
\\ Globals set by the caller before read("ptsearch.gp"):
\\   FAM, PP, QQ   family and parameter (model from ../families/fam.gp), or A = [a1,a2,a3,a4,a6] (then FAM = "curve")
\\   B             search bound (box |m| <= B s, 1 <= n <= B/s per quartic, s from the Stoll-Cremona covariant)
\\   OUT           output directory (created by the caller); files: covers.bin, q<i>.txt, hits<i>.txt, sieve<i>.log,
\\                 points.txt, progress.txt, summary.txt (i = member index)
\\   TARGET        rank target: members are searched one after another (input curve first) and the search stops
\\                 as soon as TARGET independent points are in hand (0 = search everything)
\\   MEMBERS       0 = all members (default), else a vector of member indices (searched in that order)
\\   PRUNE         1 = Cassels-kernel cosets on full-2-torsion members (default in ptsearch.sh)
\\   SQ            1 = square boxes (validation), default 0
\\   ALLCOV        1 = every coset of Sel2/delta(tors) on members without full 2-torsion (allcov.gp); KNOWN, ONLYNEW
\\ Prints one summary line.
read("../families/fam.gp"); read("pts.gp");
if(type(SQ) == "t_POL", SQ = 0);
if(type(MEMBERS) == "t_POL", MEMBERS = 0);
if(type(TARGET) == "t_POL", TARGET = 0);
if(type(ONLYNEW) == "t_POL", ONLYNEW = 0);   \\ 1: kind-3 covers only for classes outside span(torsion, KNOWN)
if(type(BSEL) == "t_POL", BSEL = 0);         \\ > 0 (with MEMBERS = 0): search members with b within factor BSEL of the best, best first
SIEVE = "search/points/bin/ptsieve";
{
  my(t0 = getwalltime(), M, CL, CV, H, S = [], allp = List(), line, tsieve = 0, tprep, order, done = List(), nhit = 0, nbad = 0);
  M = if(type(A) == "t_VEC", ellminimalmodel(ellinit(A)), minmodel(FAM, PP, QQ));
  H = heightE(M);
  CL = members(M);
  foreach(KNOWN, P, if(!ellisoncurve(M, P), error("KNOWN point not on the input curve")); listput(allp, P));   \\ seeds the independent set
  CV = [];
  order = if(MEMBERS != 0, MEMBERS, [1..#CL[2]]);
  if(MEMBERS == 0 && BSEL > 0,      \\ members with b >= max b / BSEL, largest b first (helpers/BUDGET.md)
    my(bv = vector(#CL[2], i, ellbsd(CL[2][i][1])), bm = vecmax(bv));
    for(i = 1, #bv, write(Str(OUT, "/members_b.txt"), "member ", i, " tors ", elltors(CL[2][i][1])[2], " log10b ", Strprintf("%.4f", log(bv[i]) / log(10)),
                                     " ", if(bv[i] * BSEL >= bm, "searched", "skipped")));
    order = vecsort(select(i -> bv[i] * BSEL >= bm, [1..#bv]), i -> -bv[i]));
  tprep = getwalltime() - t0;
  foreach(order, i,
    my(tc = getwalltime(), nc = #CV);
    CV = concat(CV, covers(CL, 1, [i]));            \\ covers built member by member (early stop at TARGET saves the rest)
    system(Str("rm -f ", OUT, "/covers.bin"));      \\ writebin appends; read() of several objects returns them all
    writebin(Str(OUT, "/covers.bin"), [M, CL, CV]);
    tprep += getwalltime() - tc;
    my(qs = select(q -> q > nc && CV[q][2] == i && (!ONLYNEW || CV[q][1] != 3 || CV[q][6]), [1..#CV]), ts = getwalltime(), qf = Str(OUT, "/q", i, ".txt"), hf = Str(OUT, "/hits", i, ".txt"), lf = Str(OUT, "/sieve", i, ".log"), SL, R);
    if(#qs == 0, next);
    writequartics(CV, qf, B, SQ, qs);
    system(Str(SIEVE, " ", qf, " --out ", hf, " 2> ", lf));
    tsieve += getwalltime() - ts;
    SL = readstr(lf);
    if(#SL == 0 || strsplit(SL[#SL], " ")[2] != "total", error("ptsieve failed: ", if(#SL, SL[#SL], "no log")));
    R = readhits(CL, CV, hf, qs); nbad += R[2]; nhit += #R[1] + R[2];
    foreach(R[1], h, listput(allp, h[5]);
      write(Str(OUT, "/points.txt"), "cover ", h[1], " member ", CV[h[1]][2], " kind ", CV[h[1]][1], " m ", h[2], " n ", h[3],
            " hE ", if(ellorder(M, h[5]) > 0, "tors", Strprintf("%.4f", ellheight(M, h[5]))), " P ", h[5]));
    S = indep(M, Vec(allp));
    listput(done, i);
    write(Str(OUT, "/progress.txt"), "member ", i, " covers ", #qs, " indep ", #S, " sieve_s ", (getwalltime() - ts) \ 1000);
    if(TARGET > 0 && #S >= TARGET, break));
  if(#S, S = ellsaturation(M, S, 100); S = indep(M, S));
  write(Str(OUT, "/points.txt"), "# independent ", #S);
  foreach(S, P, write(Str(OUT, "/points.txt"), "P ", P[1], " ", P[2]));
  line = Str(FAM, " ", PP, " ", QQ, " ", Strprintf("%.6f", H), " B ", B, " covers ", #select(c -> setsearch(Set(Vec(done)), c[2]), CV),
             " members ", strjoin(apply(x -> Str(x), Vec(done)), ","), " hits ", nhit, " points ", #allp, " indep ", #S,
             " heights ", strjoin(vector(#S, i, Strprintf("%.3f", ellheight(M, S[i]))), ","),
             " prep_s ", tprep \ 1000, " sieve_s ", tsieve \ 1000, " total_s ", (getwalltime() - t0) \ 1000);
  print(line); write(Str(OUT, "/summary.txt"), line);
}
quit
