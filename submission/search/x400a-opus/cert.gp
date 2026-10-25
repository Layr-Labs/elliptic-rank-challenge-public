read("lib.gp");
fmtq(x)= if(denominator(x)==1, Str(x), Str(numerator(x),"/",denominator(x)));
cert(p,q,eff,fname,hdr)={
  my(F=x400afamily(p,q), E=F[1], M, w, P, r, pts, G);
  M=ellminimalmodel(E,&w); P=apply(x->ellchangepoint(x,w),F[2]);
  r=ellrank(M,eff,P); pts=r[4]; print("rank bounds ",r[1]," ",r[2]," pts ",#pts, " H=",NH(M)," tors=",elltors(M)[2]);
  pts=ellsaturation(M,pts,100);
  G=ellheightmatrix(M,pts); print("minev ",vecmin(qfjacobi(G)[1]));
  my(s=hdr); s=Str(s,"\nE ",M.a1," ",M.a2," ",M.a3," ",M.a4," ",M.a6);
  for(i=1,#pts, s=Str(s,"\nP ",fmtq(pts[i][1])," ",fmtq(pts[i][2])));
  s=Str(s,"\n#CLAIM rank=",#pts);
  write(fname, s); #pts;
}
