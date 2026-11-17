default(realprecision,57);
read("submission/search/mestre-own/mestre.gp");
b=[-11,-10,1,5,6,9];
print("OWN_SEXTIC=",prod(i=1,6,'x-b[i]));
for(j=1,8,my(t=if(j<=5,j,if(j==6,7/3,if(j==7,11/7,13/5))),a=concat(b+vector(6,i,t),b-vector(6,i,t)),qr=mesqr(a));if(poldegree(qr[2])!=4 || poldisc(qr[2])==0,print("degenerate t=",t);next);my(F=mesPts(a,500),P=indep(F[1],F[2]));print("OWN_M12 t=",t," H=",H(F[1][1..5])," pts=",F[3]," indep=",#P));
quit;
