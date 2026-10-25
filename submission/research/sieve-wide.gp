default(realprecision,38);
ps=primes([2,251]); ws=vector(#ps,i,log(ps[i])/ps[i]); seen=Map();
{
for(f=5,9,
 N=if(f<=6,512,if(f<=8,256,160));
 for(q=1,N,
  for(p=-N,N,
   if(gcd(p,q)!=1,next);
   t=p/q;
   if(t==0 || t==1,next);
   if(f==5,b=t;c=t);
   if(f==6,c=t;b=t+t^2);
   if(f==7,c=t^2-t;b=t*c);
   if(f==8,b=(2*t-1)*(t-1);c=b/t);
   if(f==9,c=t^2*(t-1);b=c*(t^2-t+1));
   E=ellinit([1-c,-b,-b,0,0]);
   if(#E==0,next);
   E=ellminimalmodel(E);
   if(mapisdefined(seen,E.j),next);
   mapput(seen,E.j,1);
   h=log(max(abs(E.c4)^3,E.c6^2))/12;
   lim=if(f==5,8.051843,if(f==6,7.724108,if(f==7,15.121746,if(f==8,16.214667,17.847172))));
   if(h>=lim,next);
   s=0;for(i=1,#ps,if(E.disc%ps[i],s+=ellap(E,ps[i])*ws[i]));
   if(s < -8,print([f,p,q,h,s,[E.a1,E.a2,E.a3,E.a4,E.a6]]));
  );
 );
 print("FAMILY_DONE ",f," UNIQUE ",#seen);
);
}
quit;
