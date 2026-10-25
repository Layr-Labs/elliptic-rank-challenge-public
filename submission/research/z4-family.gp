\\ Elkies-Klagsbrun, published version, section 11, equations (5),(6).
z4curve(f,t)={
 my(a1,a2);
 if(f==1,a1=(8*t-1)*(32*t+7);a2=8*(t+1)*(15*t-8)*(31*t-7),
 a1=-8*(80*t+9);a2=2*(t-2)*(2*t-1)*(18*t-1)*(2*t-81));
 return(ellinit([a1,a2,a1*a2,0,0]));
};
z4family(f,t)={
 my(E=z4curve(f,t),xs,P=List(),dd,yy,xx);
 if(#E==0,return([[],[]]));
 if(f==1,xs=[-15/4*(t+1)*(31*t-7)*(32*t+7),(8*t-1)*(15*t-8)*(31*t-7)*(32*t+7),-(t+1)*(8*t-1)*(15*t-8)*(32*t+7),-4*(t+1)*(2*t+5)*(15*t-8)*(32*t+7)],
 xs=[154*(t-2)*(2*t-1)*(18*t-1),-1456*(t-2)*(2*t-1)*(2*t-81),16*(t-2)*(2*t-81)*(22*t+21),6*(2*t-5)*(t-2)*(2*t-81)*(18*t-1)]);
 for(i=1,#xs,xx=xs[i];dd=(E.a1*xx+E.a3)^2+4*(xx^3+E.a2*xx^2);if(issquare(dd,&yy),listput(P,[xx,(-E.a1*xx-E.a3+yy)/2]),error("section discriminant not square")));
 return([E,Vec(P)]);
};
