default(realprecision,57);
read("submission/search/mestre-own/mestre.gp");
b=[148,116,104,57,25,0];
for(j=1,6,my(t=if(j<=3,j,if(j==4,23550/150-150/2,if(j==5,23550/154-154/2,23550/126-126/2))),a=concat(b+vector(6,i,t),b-vector(6,i,t)),qr=mesqr(a));if(poldegree(qr[2])!=4||poldisc(qr[2])==0,next);my(F=mesPts(a,500),P=indep(F[1],F[2]));print("NAGAO t=",t," H=",H(F[1][1..5])," pts=",F[3]," indep=",#P));
quit;
