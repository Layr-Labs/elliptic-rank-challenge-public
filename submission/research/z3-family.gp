read("/home/linux/ec-rank/submission/research/pairings.gp");
z3family(rr,ss,vv,lam,kk,t) = {
 my(roots=[-vv,-ss,-rr,rr,ss,vv],Q,ks,qs,pairs,P=List(),uu,ww,xx,yy,E,a1,a3);
 a1=t^2+lam;a3=kk*(t^2-rr^2)*(t^2-ss^2)*(t^2-vv^2);
 E=ellinit([a1,0,a3,0,0]);
 if(#E==0,return([[],[]]));
 for(j=1,#matchings,
  pairs=matchings[j];
  Q=matrix(3,3,i,k,if(i==1,1,if(i==2,-roots[pairs[k][1]]-roots[pairs[k][2]],roots[pairs[k][1]]*roots[pairs[k][2]])));
  if(matdet(Q)==0,next);
  ks=Q^-1*[1,0,lam]~;
  if(prod(k=1,3,ks[k])!=kk,next);
  qs=vector(3,k,(t-roots[pairs[k][1]])*(t-roots[pairs[k][2]]));
  uu=-ks[1]*qs[1];ww=ks[2]*qs[2];
  xx=uu*ww;yy=uu^2*ww;
  if(ellisoncurve(E,[xx,yy]),listput(P,[xx,yy]),error("factor-section identity failed"));
 );
 return([E,Vec(Set(P))]);
};
