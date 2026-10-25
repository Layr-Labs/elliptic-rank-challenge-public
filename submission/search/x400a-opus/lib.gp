read("x400a-family.gp");
NH(E) = my(c4=E.c4,c6=E.c6); log(max(abs(c4)^3, c6^2)*1.)/12;
\\ class members: returns list of [H, torsion, minimal curve] over the 2-power isogeny class (via ellisomat)
cls(E)={ my(M=ellisomat(E,2), L=M[1], R=List()); for(i=1,#L, my(C=ellminimalmodel(ellinit(L[i][1]))); listput(R,[NH(C), elltors(C)[2], C])); Vec(R)};
\\ L=[[p,q],...]: min H over class members with torsion [2,2]; rank on E
st(L, Hmax, eff, out)={
  for(i=1,#L, my(p=L[i][1], q=L[i][2], F, E, M, w, P, cl, hb=99, tb, h0, r, t0=getabstime());
    F=iferr(x400afamily(p,q),e,[[],[]]); E=F[1]; if(#E==0,next);
    M=ellminimalmodel(E,&w); h0=NH(M);
    cl=cls(M); for(j=1,#cl, if(cl[j][2]==[2,2] && cl[j][1]<hb, hb=cl[j][1]));
    if(hb>=Hmax, next);
    P=apply(x->ellchangepoint(x,w),F[2]);
    r=iferr(ellrank(M,eff,P),e,[-1,-1,0,[]]);
    write(out, strprintf("%d %d %.6f %.6f %d %d %d %d %d", p, q, hb, h0, r[1], r[2], #r[4], #cl, getabstime()-t0)));
}
