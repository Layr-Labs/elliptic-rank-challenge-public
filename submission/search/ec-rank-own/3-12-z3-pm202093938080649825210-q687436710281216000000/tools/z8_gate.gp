\\ Owner late stage for z8 candidates. Variables set by the driver: infile, target, outfile.
default(parisizemax, 4000000000);
read("work/rank/rankclass.gp");
e0(p,q) = { my(u=2*p-q, w=p*(p-q)); ellinit([0, u^4-8*w^2, 0, 16*w^4, 0]); }
e1(p,q) = { my(u=2*p-q, w=p*(p-q), Q=u^2+4*w); ellinit([0, -(u^4+q^2*Q), 0, u^4*q^2*Q, 0]); }
{
  my(L = readstr(infile), n = 0, pass = 0, t0 = getabstime());
  for(i = 1, #L,
    my(v = strsplit(L[i], " "), p, q, E1, r, E0, H, rc);
    if(#v < 3, next);
    p = eval(v[2]); q = eval(v[3]); n++;
    E1 = ellminimalmodel(e1(p,q));
    r = ellrank(E1);
    if(r[2] < target, next);
    pass++;
    E0 = ellminimalmodel(e0(p,q));
    H = log(max(abs(E0.c4)^3, E0.c6^2))/12;
    rc = rankclass(E0, 1, target);
    write(outfile, "z8 ", p, " ", q, " H ", precision(H,9), " r2E1 ", r[2], " lo ", rc[1], " hi ", rc[2], " npts ", #rc[3], " ms ", rc[4], " ainv ", E0[1..5], " pts ", rc[3]);
  );
  write(outfile, "# done ", infile, " n=", n, " pass_r2>=", target, ": ", pass, " ms=", getabstime()-t0);
}
quit
