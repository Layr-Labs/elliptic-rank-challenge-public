{
forprime(pr=2,2000,
 v=vector(pr);
 for(r=0,pr-1,t=Mod(r,pr);
  A=37-84*t+102*t^2-36*t^3-3*t^4;
  B=32*(t-1)^3*(t+1)^3*(3*t-5);
  E=ellinit([0,A,0,B,0]);v[r+1]=if(#E,pr+1-ellcard(E),9999);
 );
 print(pr," ",v);
);
}
quit;
