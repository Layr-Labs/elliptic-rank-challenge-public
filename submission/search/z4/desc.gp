\\ 2-isogeny descent point search on E (needs rational 2-torsion); returns points on E
sqfdivs(B)={my(f=factor(abs(B))[,1]~, L=List()); forsubset(#f,s, my(d=prod(i=1,#s,f[s[i]])); listput(L,d); listput(L,-d)); Vec(L)};
isogside(A,B,bnd,maxd)={my(pts=List(),D=sqfdivs(B));
  if(#D>maxd, D=D[1..maxd]);
  foreach(D,d, my(q=d*'x^4+A*'x^2+B/d); if(poldisc(q)==0,next);
    my(R=iferr(hyperellratpoints(q,bnd),e,[]));
    foreach(R,r, my(M=r[1],w=r[2]); if(M==0,next); listput(pts,[d*M^2, d*M*w])));
  Vec(pts)};
desc2(E,bnd,maxd=100000)={my(v1=[1/2,0,-E.a1/2,-E.a3/2]); my(E2=ellchangecurve(E,v1));
  my(rt=nfroots(,'x^3+E2.a2*'x^2+E2.a4*'x+E2.a6)); if(#rt==0,return([]));
  my(v2=[1,rt[1],0,0]); my(E3=ellchangecurve(E2,v2)); my(A=E3.a2,B=E3.a4);
  my(P1=isogside(A,B,bnd,maxd)); my(Ap=-2*A,Bp=A^2-4*B); my(P2=isogside(Ap,Bp,bnd,maxd));
  my(out=List());
  my(back(P)=ellchangepointinv(ellchangepointinv(P,v2),v1));
  foreach(P1,P, if(ellisoncurve(E3,P), listput(out,back(P))));
  foreach(P2,P, my(x=P[1],y=P[2]); my(Q=[y^2/(4*x^2), y*(x^2-Bp)/(8*x^2)]); if(ellisoncurve(E3,Q), listput(out,back(Q))));
  Vec(out)};
indep(E,pts)={my(S=List()); foreach(pts,P, if(!ellisoncurve(E,P)||ellorder(E,P)>0,next); my(c=concat(Vec(S),[P])); my(M=ellheightmatrix(E,c)); my(ev=vecmin(real(qfjacobi(M)[1]))); if(ev>1e-6, listput(S,P))); Vec(S)};
