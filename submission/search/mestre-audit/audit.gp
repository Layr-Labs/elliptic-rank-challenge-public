default(realprecision,57);
read("submission/search/mestre-own/mestre.gp");
qrpoly(p,n)={my(s=Ser(subst(p,'x,1/'x)*'x^(2*n),'x),sq,q); sq=sqrt(s); q=sum(i=0,n,polcoef(sq,i,'x)*'x^(n-i)); [q,q^2-p]};
F='x^6+'a*'x^4+'b*'x^3+'c*'x^2+'d*'x+'e;
z=qrpoly(subst(F,'x,'x-'t)*subst(F,'x,'x+'t),6);
print("TWELVE_ROOT_R5=",polcoef(z[2],5,'x));
print("TWELVE_ROOT_R4=",polcoef(z[2],4,'x));
for(j=1,4,my(a=if(j==1,[0,1,7,8,15,16,24,25,34,35],if(j==2,[0,2,8,10,20,22,30,32,40,42],if(j==3,[0,3,12,14,26,27,39,42,55,57],[0,4,8,12,20,28,36,40,48,53]))),qr=mesqr(a),r=qr[2],A=polcoef(r,4),B=polcoef(r,3),C=polcoef(r,2),D=polcoef(r,1),E=polcoef(r,0),I=12*A*E-3*B*D+C^2,J=72*A*C*E+9*B*C*D-27*A*D^2-27*E*B^2-2*C^3); if(4*I^3-J^2==0,next);my(M=mesPts(a,200));my(P=indep(M[1],M[2]));print("TEN_TEST=",a," H=",H(M[1][1..5])," quartic_pts=",M[3]," indep=",#P," qden=",denominator(qr[1])));
quit;
