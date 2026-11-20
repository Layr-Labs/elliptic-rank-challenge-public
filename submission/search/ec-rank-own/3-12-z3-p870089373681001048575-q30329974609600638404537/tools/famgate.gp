\\ Owner generic late stage. Driver sets: infile, target, outfile.  Input lines "fam p q ...".
\\ For each line whose (p,q) is the canonical representative of its symmetry orbit: minimal model, ellrank (effort 0);
\\ if the upper bound r2 >= target: rankclass over the isogeny class (effort 1) and write a result line.
default(parisizemax, 4000000000);
read("search/families/fam.gp");
read("work/rank/rankclass.gp");
{
  my(L = readstr(infile), n = 0, ncanon = 0, pass = 0, t0 = getabstime(), hist = vector(12));
  for(i = 1, #L,
    my(v = strsplit(L[i], " "), f, p, q, E, r, H, rc, c);
    if(#v < 3, next);
    f = v[1]; p = eval(v[2]); q = eval(v[3]); n++;
    c = canon(f, p, q);
    if(c != [p, q], next);
    ncanon++;
    E = minmodel(f, p, q);
    r = ellrank(E);
    hist[min(r[2], 11) + 1]++;
    if(r[2] < target, next);
    pass++;
    H = heightE(E);
    rc = rankclass(E, 1, target);
    write(outfile, f, " ", p, " ", q, " H ", precision(H, 9), " r2 ", r[2], " lo ", rc[1], " hi ", rc[2], " npts ", #rc[3], " ms ", rc[4], " ainv ", E[1..5], " pts ", rc[3]);
  );
  write(outfile, "# done ", infile, " lines=", n, " canonical=", ncanon, " r2>=", target, ": ", pass, " r2hist(0..11)=", hist, " ms=", getabstime() - t0);
}
quit
