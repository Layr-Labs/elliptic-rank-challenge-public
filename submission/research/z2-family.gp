\\ Elkies-Klagsbrun, arXiv:2003.00077, equations (3),(4),(7).
\\ Family formulas only; all specializations are chosen by this search.
z2family(u,t) = {
 my(A,B,b1,b2,b3,b4,b5,b6,b7,b8,xs,E,P=List(),yy,x);
 A=(u^8-18*u^6+163*u^4-1152*u^2+4096)*t^4+(3*u^7-35*u^5-120*u^3+1536*u)*t^3+(u^8-13*u^6+32*u^4-152*u^2+1536)*t^2+(u^7+3*u^5-156*u^3+672*u)*t+(3*u^6-33*u^4+112*u^2-80);
 b1=(u^2+u-8)*t-u+2;
 b2=(u^2-u-8)*t-u-2;
 b3=(u^2-u-8)*t+u^2+u-10;
 b4=(u^2+u-8)*t-u^2+u+10;
 b5=(u^2-7*u+8)*t-u^2+u+2;
 b6=(u^2+7*u+8)*t+u^2+u-2;
 b7=(u^2+5*u+8)*t+u^2+3*u+2;
 b8=(u^2-5*u+8)*t-u^2+3*u-2;
 B=b1*b2*b3*b4*b5*b6*b7*b8;
 E=ellinit([0,2*A,0,B,0]);
 if(#E==0,return([[],[]]));
 xs=[-b1*b2*b3*b6,-b1*b2*b4*b5,4*b1*b2*b5*b6,b1*b3*b4*b6,-b1*b3*b4*b7,b1*b3*b4*b8,b1*b3*b5*b6,-b1*b5*b6*b7];
 for(i=1,#xs,x=xs[i];if(issquare(x^3+2*A*x^2+B*x,&yy),listput(P,[x,yy])));
 return([E,Vec(P)]);
};
