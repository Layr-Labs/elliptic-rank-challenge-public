{
forprime(p=2,1000,
 for(f=5,6,
  v=vector(p);
  for(r=0,p-1,
   c=r;b=if(f==5,r,r+r^2);
   E=ellinit([1-c,-b,-b,0,0]);
   v[r+1]=if(#E && E.disc%p,ellap(E,p),0);
  );
  print(f," ",p," ",v);
 );
);
}
quit;
