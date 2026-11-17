u='u;t='t;
AA=(u^8-18*u^6+163*u^4-1152*u^2+4096)*t^4+(3*u^7-35*u^5-120*u^3+1536*u)*t^3+(u^8-13*u^6+32*u^4-152*u^2+1536)*t^2+(u^7+3*u^5-156*u^3+672*u)*t+3*u^6-33*u^4+112*u^2-80;
ff=[(u^2+u-8)*t-u+2,(u^2-u-8)*t+u^2+u-10,(u^2-7*u+8)*t-u^2+u+2,(u^2+5*u+8)*t+u^2+3*u+2];
BB=vector(8,k,if(k%2,ff[(k+1)/2],-subst(subst(ff[k/2],t,-t),u,-u)));
A=2*AA;B=prod(k=1,8,BB[k]);
XX=[-BB[1]*BB[2]*BB[3]*BB[6],-BB[1]*BB[2]*BB[4]*BB[5],4*BB[1]*BB[2]*BB[5]*BB[6],BB[1]*BB[3]*BB[4]*BB[6],-BB[1]*BB[3]*BB[4]*BB[7],BB[1]*BB[3]*BB[4]*BB[8],BB[1]*BB[3]*BB[5]*BB[6],-BB[1]*BB[5]*BB[6]*BB[7]];
{
for(j=1,2,my(m=if(j==1,2,-3),uv=2*(m^2-m-1)/(m^2+1),av=subst(A,u,uv),bv=subst(B,u,uv),xs=concat(subst(XX,u,uv),[-(m-1)^2*subst(BB[1]*BB[2]*BB[3]*BB[8],u,uv)]));for(k=1,9,if(!issquare(xs[k]^3+av*xs[k]^2+bv*xs[k]),error("section wrong ",j," ",k)));print("SECTIONS=",[uv,#xs]);print("AC=",vector(5,k,polcoef(av*625^2,k-1)));print("BC=",vector(9,k,polcoef(bv*625^4,k-1)));
for(i=1,4,my(tt=[1/2,-1/3,2/3,-3/5][i],E=ellinit([0,subst(av,t,tt),0,subst(bv,t,tt),0]),M=ellminimalmodel(E,&ch),pts=vector(9,k,my(xx=subst(xs[k],t,tt));ellchangepoint([xx,ellordinate(E,xx)[1]],ch)),R=ellrank(M,1,pts));print("TEST=",[uv,tt,log(max(abs(M.c4)^3,M.c6^2))/12.,R[1],R[2],#R[4]])));
}
print("COMPLETE");quit;
