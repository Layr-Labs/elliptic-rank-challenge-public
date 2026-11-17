default(realprecision,57);
read("submission/search/mestre-own/mestre.gp");
b=[-17,-16,10,11,14,17];
for(j=1,6,my(t=if(j<=4,j,if(j==5,1/2,3/2)),a=concat(b+vector(6,i,t),b-vector(6,i,t)),qr=mesqr(a));print("M12_T=",t," deg=",poldegree(qr[2]));if(poldegree(qr[2])!=4,next);my(F=mesPts(a,3000),P=indep(F[1],F[2]));print("M12_RESULT t=",t," H=",H(F[1][1..5])," pts=",F[3]," indep=",#P));
quit;
