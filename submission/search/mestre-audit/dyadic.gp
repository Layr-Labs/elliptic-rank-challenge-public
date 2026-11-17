default(realprecision,38);
read("submission/search/mestre-own/mestre.gp");
a=[0,11,17,33,47,51,59,62,65,68];q0=mesqr(a);q4=mesqr(4*a);F=mesPts(a,100);G=mesPts(4*a,100);print("DYADIC denominator=",lcm(apply(denominator,Vec(q0[1])))," after4=",lcm(apply(denominator,Vec(q4[1])))," same_j=",F[1].j==G[1].j," same_model=",F[1][1..5]==G[1][1..5]," H=",H(F[1][1..5]));print("SCALE_Q_IDENTITY=",q4[1]==4^5*subst(q0[1],x,x/4));print("SCALE_R_IDENTITY=",q4[2]==4^10*subst(q0[2],x,x/4));
quit;
