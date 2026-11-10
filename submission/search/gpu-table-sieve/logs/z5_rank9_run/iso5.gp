\\ 5-isogeny partner of the Z/5 board curve via Velu on <P5> (P5 a point of order 5)
default(parisizemax, 2*10^9);
Hn(E)=log(max(abs(E.c4)^3, E.c6^2))/12.;
minim(a)={
  my(E=ellinit(a), g=gcd(E.c4,E.c6), F, lr);
  if(g==1, return(E));
  F=factor(g)[,1];
  for(i=1,#F, lr=elllocalred(E,F[i]); E=ellinit(ellchangecurve(E,lr[3])[1..5]));
  E
};
E = ellinit([1,0,0,-6612038604938639760365630912643765, 207009790989398164552588533643664645493811922402721]);
T = elltors(E); P5 = T[3][1];
print("E: torsion ", T[2], " H=", Hn(E));
iso = ellisogeny(E, P5);
Q = minim(iso[1][1..5]);
print("E/<P5>: H=", Hn(Q), " torsion=", elltors(Q)[2], " c4=", Q.c4, " c6=", Q.c6);
\\ the quotient's own isogenies: dual returns E; other rational 5-isogenies of E? check the 5-division polynomial factor of degree 2 (kernel polynomials of degree 2)
D5 = elldivpol(E,5); print("degree of 5-division polynomial: ", poldegree(D5)); F5 = factor(D5)[,1]; print("degrees of its irreducible factors over Q: ", vector(#F5,i,poldegree(F5[i])));
