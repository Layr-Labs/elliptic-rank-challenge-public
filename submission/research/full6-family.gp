\\ Dujella-Peral, arXiv:1306.0027, universal Z/2 x Z/6 model.
full6curve(p,q)={
 my(A=37*q^4-84*p*q^3+102*p^2*q^2-36*p^3*q-3*p^4,
 B=32*q*(p-q)^3*(p+q)^3*(3*p-5*q));
 return(ellinit([0,A,0,B,0]));
};
