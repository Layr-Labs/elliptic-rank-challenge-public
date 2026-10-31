\\ Mestre ten-root construction: helpers used for the 1/19 search (PARI/GP 2.17.4). Standalone.
H(e)=my(E=ellinit(e));log(max(abs(E.c4)^3,E.c6^2))/12;
readsplit(fn)={my(V=readstr(fn)); vector(#V,i,eval(strsplit(V[i]," ")))};
\\ p = prod(x-a_i) = q^2 - r with q monic quintic, deg r <= 4
mesqr(a)={my(p=prod(i=1,#a,'x-a[i]),n=#a/2,q);
 my(s=Ser(subst(p,'x,1/'x)*'x^(2*n),'x)); my(sq=sqrt(s)); q=sum(i=0,n,polcoef(sq,i,'x)*'x^(n-i)); [q, q^2-p]};
\\ rational points on the quartic y^2=r(x) (hyperellratpoints), mapped to the Weierstrass model
\\ by Washington, Elliptic Curves (2nd ed.), Thm 2.17, with base point (a_1, q(a_1)); returns [minimal model, points, #quartic points]
mesPts(a,bnd)={my(qr=mesqr(a),Q=qr[1],r=qr[2],x0=a[1],y0=subst(Q,'x,x0));
 my(R=subst(r,'x,'x+x0)); my(A=polcoef(R,4),B=polcoef(R,3),C=polcoef(R,2),D=polcoef(R,1),q=y0);
 my(a1=D/q,a2=C-D^2/(4*q^2),a3=2*q*B,a4=-4*q^2*A,a6=a2*a4); my(E=ellinit([a1,a2,a3,a4,a6]));
 my(S=hyperellratpoints(r,bnd)); my(P=List());
 foreach(S,pt, my(u=pt[1]-x0,v=pt[2]); if(u==0,next);
   my(X=(2*q*(v+q)+D*u)/u^2, Y=(4*q^2*(v+q)+2*q*(D*u+C*u^2)-D^2*u^2/(2*q))/u^3);
   if(ellisoncurve(E,[X,Y]), listput(P,[X,Y])));
 my(w,M=ellinit(ellminimalmodel(E,&w))); [M, apply(z->ellchangepoint(z,w),Vec(P)), #S]};
indep(M,P)={my(S=List()); P=vecsort(P,z->ellheight(M,z)); for(i=1,#P, my(T=concat(Vec(S),[P[i]])); if(matdet(ellheightmatrix(M,T))>1e-8*10^(#T), listput(S,P[i]))); Vec(S)};
\\ quick rank lower bound for sieve survivors: file lines "a1..a10 H S"
quickMes(fn,bnd)={my(L=readsplit(fn)); for(i=1,#L, my(a=L[i][1..10]); my(F=iferr(mesPts(a,bnd),e,0)); if(F==0,next); my(M=F[1]); my(I=indep(M,F[2]));
  printf("%s H=%.6f S=%.3f r>=%d\n", strjoin(apply(z->Str(z),a)," "), H(M[1..5]), L[i][#L[i]], #I))};
fmtq(x)= if(denominator(x)==1, Str(x), Str(numerator(x),"/",denominator(x)));
certMes(a,r,fn,hdr)={my(F=mesPts(a,10^5),M=F[1],P=F[2]);
 P=ellsaturation(M,indep(M,P),100); my(S=indep(M,P)); if(#S<r, print("only ",#S); return(0)); S=S[1..r];
 my(G=ellheightmatrix(M,S)); print("H=",H(M[1..5])," minev=",vecmin(qfjacobi(G)[1])," tors=",elltors(M)[2]);
 my(s=hdr); s=Str(s,"\nE ",M.a1," ",M.a2," ",M.a3," ",M.a4," ",M.a6);
 for(i=1,#S, s=Str(s,"\nP ",fmtq(S[i][1])," ",fmtq(S[i][2]))); s=Str(s,"\n#CLAIM rank=",r);
 my(f=fileopen(fn,"w")); filewrite(f,s); fileclose(f); 1};
