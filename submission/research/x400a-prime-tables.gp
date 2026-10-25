{
forprime(pr=2,2000,
 v=vector(pr);
 for(r=0,pr-1,t=Mod(r,pr);A=(121-99*t)*(26*t-23)*(75*t-44);B=t*(21-74*t)*(25*t-36)*(t+77);E=ellinit([0,A+B,0,A*B,0]);v[r+1]=if(#E,pr+1-ellcard(E),9999));
 print(pr," ",v)
);
}
quit;
