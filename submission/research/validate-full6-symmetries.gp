A=37-84*x+102*x^2-36*x^3-3*x^4;
B=32*(x-1)^3*(x+1)^3*(3*x-5);
C4=16*(A^2-3*B);C6=-64*A^3+288*A*B;
if(subst(C4,x,2-x)!=C4 || subst(C6,x,2-x)!=C6,error("reflection invariants"));
M=(7-x)/(5-3*x);um=16/(5-3*x)^2;
if(numerator(subst(C4,x,M)-C4*um^4)!=0 || numerator(subst(C6,x,M)-C6*um^6)!=0,error("order-three invariants"));
print("FULL6_SYMMETRIES_EXACT");
quit;
