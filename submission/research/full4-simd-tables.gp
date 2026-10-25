{
forprime(p=2,2000,
 v=vector(p);
 for(r=0,p-1,t=Mod(r,p);E=ellinit([0,2*(t^2+1),0,(t^2-1)^2,0]);v[r+1]=if(#E,p+1-ellcard(E),9999));
 print(p," ",v);
);
}
quit;
