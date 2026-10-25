read("submission/research/full6-family.gp");
{
for(r=0,127,
 for(s=0,127,
  if(r%2==0 && s%2==0,next);
  pp=r+512;qq=s+256;E=full6curve(pp,qq);
  while(#E==0,pp+=128;E=full6curve(pp,qq));
  red=elllocalred(E,2);
  print(r," ",s," ",valuation(red[3][1],2)," ",red[1]);
 );
);
}
quit;
