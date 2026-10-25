read("submission/research/z4-family.gp");
{
for(f=1,2,
 forprime(p=2,1000,vap=vector(p);
  for(j=0,p-1,t=Mod(j,p);E=z4curve(f,t);vap[j+1]=if(#E,p+1-ellcard(E),9999));
  print(f," ",p," ",concat(vap,[9999]));
 );
);
}
quit;
