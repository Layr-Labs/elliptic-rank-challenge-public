read("submission/research/z2-family.gp");
default(realprecision,38);
ps=primes([2,1000]);ws=vector(#ps,i,log(ps[i])/ps[i]);seen=Map();
{
for(v=1,20,
 for(w=1,16,
  if(gcd(v,w)!=1,next);
  u=v/w;
  if(u==1 || u==2,next);
  for(q=1,30,
   for(p=-30,30,
    if(gcd(p,q)!=1,next);
    F=z2family(u,p/q);
    if(#F[1]==0 || #F[2]!=8,next);
    origE=F[1];
    E=ellminimalmodel(origE);
    E2=ellminimalmodel(ellinit([0,-2*origE.a2,0,origE.a2^2-4*origE.a4,0]));
    h=log(max(abs(E.c4)^3,E.c6^2))/12;
    h2=log(max(abs(E2.c4)^3,E2.c6^2))/12;
    iso=0;if(h2<h,h=h2;E=E2;iso=1);
    if(h>=14.5,next);
    if(mapisdefined(seen,E.j),next);
    mapput(seen,E.j,1);
    s=0;for(i=1,#ps,if(E.disc%ps[i],s+=ellap(E,ps[i])*ws[i]));
    print([2,p,q,h,s,[E.a1,E.a2,E.a3,E.a4,E.a6],Str(u),iso]);
   );
  );
 );
 print("U_DONE ",v," RETAINED ",#seen);
);
}
quit;
