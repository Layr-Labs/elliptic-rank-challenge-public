default(realprecision,38);
ps=primes([2,1000]); ws=vector(#ps,i,log(ps[i])/ps[i]); seen=Map();
{
forstep(f=10,12,2,
 for(q=1,96,
  for(p=-96,96,
   if(gcd(p,q)!=1,next);
   t=p/q;
   if(t==0 || t==1,next);
   if(f==10,c=(2*t^3-3*t^2+t)/(t-(t-1)^2);b=c*t^2/(t-(t-1)^2));
   if(f==12,c=(3*t^2-3*t+1)*(t-2*t^2)/(t-1)^3;b=c*(2*t-2*t^2-1)/(t-1));
   E=ellinit([1-c,-b,-b,0,0]);
   if(#E==0,next);
   E=ellminimalmodel(E);
   if(mapisdefined(seen,E.j),next);
   mapput(seen,E.j,1);
   h=log(max(abs(E.c4)^3,E.c6^2))/12;
   lim=if(f==10,20.305266,27.243549);
   if(h>=lim,next);
   s=0;for(i=1,#ps,if(E.disc%ps[i],s+=ellap(E,ps[i])*ws[i]));
   print([f,p,q,h,s,[E.a1,E.a2,E.a3,E.a4,E.a6]]);
  );
 );
);
}
quit;
