read("submission/research/z2-family.gp");
default(realprecision,38);
ps=primes([2,1000]);ws=vector(#ps,i,log(ps[i])/ps[i]);seen=Map();
{
for(v=1,7,
 for(w=1,8,
  if(gcd(v,w)!=1,next);
  u=v/w;
  if(u==1 || u==2 || u==4/7,next);
  for(q=1,20,
   for(p=-20,20,
    if(gcd(p,q)!=1,next);
    F=z2family(u,p/q);
    if(#F[1]==0 || #F[2]!=8,next);
    E=ellminimalmodel(F[1]);
    if(mapisdefined(seen,E.j),next);
    mapput(seen,E.j,1);
    h=log(max(abs(E.c4)^3,E.c6^2))/12;
    if(h>=15.2,next);
    s=0;for(i=1,#ps,if(E.disc%ps[i],s+=ellap(E,ps[i])*ws[i]));
    print([2,p,q,h,s,[E.a1,E.a2,E.a3,E.a4,E.a6],Str(u)]);
   );
  );
 );
);
}
quit;
