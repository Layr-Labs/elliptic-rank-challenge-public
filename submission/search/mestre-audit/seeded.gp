default(realprecision,57);
read("submission/search/mestre-own/mestre.gp");
\\ Previous bounds omit rational seed points; append explicitly.
seedmes(a,extra=[])={my(qr=mesqr(a),Q=qr[1],r=qr[2],x0=a[1],q=subst(Q,'x,x0),R=subst(r,'x,'x+x0),A=polcoef(R,4),B=polcoef(R,3),C=polcoef(R,2),D=polcoef(R,1),E=ellinit([D/q,C-D^2/(4*q^2),2*q*B,-4*q^2*A,(C-D^2/(4*q^2))*(-4*q^2*A)]),w,M=ellinit(ellminimalmodel(E,&w)),P=List(),S=concat(concat(apply(xx->[xx,subst(Q,'x,xx)],a),apply(xx->[xx,-subst(Q,'x,xx)],a)),extra));listput(P,ellchangepoint([-E.a2,E.a1*E.a2-E.a3],w));my(sqA);if(issquare(A,&sqA),listput(P,ellchangepoint([2*q*sqA,0],w));listput(P,ellchangepoint([-2*q*sqA,0],w)));foreach(S,pt,my(u=pt[1]-x0,v=pt[2]);if(u==0,next);my(X=(2*q*(v+q)+D*u)/u^2,Y=(4*q^2*(v+q)+2*q*(D*u+C*u^2)-D^2*u^2/(2*q))/u^3);if(!ellisoncurve(E,[X,Y]),print("BAD pt=",pt," q=",q," quartic_residual=",v^2-subst(r,'x,pt[1]));error("bad map"));listput(P,ellchangepoint([X,Y],w)));[M,Vec(P)]};
b=[148,116,104,57,25,0];
for(j=1,3,my(t=if(j==1,7/2,if(j==2,239/2,-83/154)),a=concat(b+vector(6,i,t),b-vector(6,i,t)),extra=[[(t+703)/15,1248*t*(-224*t^3-844*t^2+900484*t+2161725)/75]],F=seedmes(a,extra),P=indep(F[1],F[2]));print("NAGAO_SEEDED t=",t," H=",H(F[1][1..5])," count=",#P," det=",matdet(ellheightmatrix(F[1],P))));
quit;
