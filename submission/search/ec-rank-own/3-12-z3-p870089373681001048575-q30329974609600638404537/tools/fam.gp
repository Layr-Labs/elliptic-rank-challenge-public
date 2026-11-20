\\ fam.gp -- Kubert/Tate torsion families over Q, homogeneous integral models in (p,q), gcd(p,q)=1.
\\ Families: z7 z8 z9 z10 z12 z2z6 z2z8.  Tate normal form E(b,c): y^2+(1-c)xy-by = x^3-bx^2, P=(0,0), d=p/q.
\\ z2z8 uses parameter a=p/q with d = a(8a+2)/(8a^2-1) in the Z/8 formulas.
\\ Models derived by derive2.gp (scale u chosen place by place on P^1 to make c4,c6,disc integral and of
\\ minimal valuation at every irreducible form); z2z8 additionally scaled by lambda=2 to clear a 1/2.
\\ API:  FAMS, famidx(f), torsN(f), torsG(f), curve(f,p,q), minmodel(f,p,q), height(f,p,q), heightE(E),
\\       forms(f), syms(f), symorbit(f,p,q), canon(f,p,q), tbc(f,x), params(f,E), discconst(f)
FAMS = ["z7","z8","z9","z10","z12","z2z6","z2z8"];
famidx(f) = { for(i=1,#FAMS, if(FAMS[i]==f, return(i))); error("unknown family ", f) };
\\ elltors()[2] expected
torsG(f) = [[7],[8],[9],[10],[12],[6,2],[8,2]][famidx(f)];
torsN(f) = torsG(f)[1];

tnf(b,c) = [1-c, -b, -b, 0, 0];
\\ (b,c) of the Tate normal form as rational functions of the family parameter x
tbc(f,x) = {
  my(i=famidx(f), b, c);
  if(i==1, b=x^3-x^2; c=x^2-x,
     i==2, b=(2*x-1)*(x-1); c=b/x,
     i==3, c=x^2*(x-1); b=c*(x*(x-1)+1),
     i==4, c=(2*x^3-3*x^2+x)/(x-(x-1)^2); b=c*x^2/(x-(x-1)^2),
     i==5, c=(3*x^2-3*x+1)*(x-2*x^2)/(x-1)^3; b=c*(2*x-2*x^2-1)/(x-1),
     i==6, c=(10-2*x)/(x^2-9); b=c+c^2,
     i==7, my(e=x*(8*x+2)/(8*x^2-1)); b=(2*e-1)*(e-1); c=b/e);
  [b,c]
};

\\ homogeneous integral models (a1,a2,a3,a4,a6), see derive2.gp output
{MODELS = [
 [-p^2 + q*p + q^2, -q*p^3 + q^2*p^2, -q^3*p^3 + q^4*p^2, 0, 0],
 [-2*p^2 + 4*q*p - q^2, -2*p^4 + 3*q*p^3 - q^2*p^2, -2*q*p^5 + 3*q^2*p^4 - q^3*p^3, 0, 0],
 [-p^3 + q*p^2 + q^3, -q*p^5 + 2*q^2*p^4 - 2*q^3*p^3 + q^4*p^2, -q^4*p^5 + 2*q^5*p^4 - 2*q^6*p^3 + q^7*p^2, 0, 0],
 [2*p^3 - 2*q*p^2 - 2*q^2*p + q^3, -2*q*p^5 + 3*q^2*p^4 - q^3*p^3, -2*q^2*p^7 + 9*q^3*p^6 - 12*q^4*p^5 + 6*q^5*p^4 - q^6*p^3, 0, 0],
 [6*p^4 - 8*q*p^3 + 2*q^2*p^2 + 2*q^3*p - q^4, -12*p^8 + 54*q*p^7 - 106*q^2*p^6 + 119*q^3*p^5 - 83*q^4*p^4 + 36*q^5*p^3 - 9*q^6*p^2 + q^7*p, -12*q*p^11 + 90*q^2*p^10 - 304*q^3*p^9 + 611*q^4*p^8 - 812*q^5*p^7 + 748*q^6*p^6 - 485*q^7*p^5 + 219*q^8*p^4 - 66*q^9*p^3 + 12*q^10*p^2 - q^11*p, 0, 0],
 [p^2 + 2*q*p - 19*q^2, 2*q*p^3 - 14*q^2*p^2 + 22*q^3*p - 10*q^4, 2*q*p^5 - 14*q^2*p^4 + 4*q^3*p^3 + 116*q^4*p^2 - 198*q^5*p + 90*q^6, 0, 0],
 [2*(32*p^4 - 12*q^2*p^2 - 4*q^3*p - 1/2*q^4), 4*(-256*q*p^7 - 384*q^2*p^6 - 240*q^3*p^5 - 80*q^4*p^4 - 14*q^5*p^3 - q^6*p^2), 8*(-8192*q*p^11 - 14336*q^2*p^10 - 9728*q^3*p^9 - 2688*q^4*p^8 + 256*q^5*p^7 + 416*q^6*p^6 + 128*q^7*p^5 + 18*q^8*p^4 + q^9*p^3), 0, 0]
];}
\\ cusp forms [form(p,q), width = exponent in the minimal discriminant at primes l >= 5 off the bad set]
{FORMS = [
 [[q,7],[p,7],[p-q,7],[p^3-8*p^2*q+5*p*q^2+q^3,1]],
 [[p,8],[p-q,8],[2*p-q,4],[q,2],[8*p^2-8*p*q+q^2,1]],
 [[p,9],[q,9],[p-q,9],[p^2-p*q+q^2,3],[p^3-6*p^2*q+3*p*q^2+q^3,1]],
 [[p,10],[p-q,10],[2*p-q,5],[q,5],[p^2-3*p*q+q^2,2],[4*p^2-2*p*q-q^2,1]],
 [[p,12],[p-q,12],[2*p-q,6],[q,2],[2*p^2-2*p*q+q^2,3],[3*p^2-3*p*q+q^2,4],[6*p^2-6*p*q+q^2,1]],
 [[p-q,6],[p-5*q,6],[q,6],[p-9*q,2],[p-3*q,2],[p+3*q,2]],
 [[p,8],[q,8],[2*p+q,8],[4*p+q,8],[8*p^2+4*p*q+q^2,4],[8*p^2-q^2,2],[8*p^2+8*p*q+q^2,2]]
];}
forms(f) = FORMS[famidx(f)];
model(f,P,Q) = substvec(MODELS[famidx(f)], ['p,'q], [P,Q]);
curve(f,P,Q) = ellinit(model(f,P,Q));

\\ minimal model without factoring the discriminant: a prime can be non-minimal (or additive) only if it divides
\\ gcd(c4,c6) (or is 2,3); such primes divide fixed resultants, so gcd(c4,c6) factors instantly.
\\ The returned model is normalised like ellminimalmodel (a1,a3 in {0,1}, a2 in {-1,0,1}).
normalize(E) = {
  my(a1=E.a1, a2=E.a2, a3=E.a3, s, r, t);
  s = -(a1 - (a1%2))/2; a1 = a1+2*s;
  r = -round((a2 - s*E.a1 - s^2)/3); \\ a2' = a2 - s*a1 + 3r - s^2
  my(a2n = E.a2 - s*E.a1 + 3*r - s^2);
  if(a2n > 1, r -= 1); if(a2n < -1, r += 1);
  my(a3n = E.a3 + r*E.a1);
  t = -(a3n - (a3n%2))/2;
  ellchangecurve(E, [1, r, s, t]);
};
minE(E) = {
  my(g = gcd(E.c4, E.c6), P = Set(concat([2,3], factor(g)[,1]~)));
  foreach(P, l, my(v = elllocalred(E, l)[3]); if(v != [1,0,0,0], E = ellchangecurve(E, v)));
  normalize(E)
};
minmodel(f,P,Q) = minE(curve(f,P,Q));
heightE(E) = log(max(abs(E.c4)^3, E.c6^2))/12;
height(f,P,Q) = heightE(minmodel(f,P,Q));

\\ ---- Tate normal form of (E, P) and family parameters of a curve ----
tate(E,P) = {
  my(E2 = ellchangecurve(E,[1,P[1],0,P[2]]), s, E3, u, E4);
  s = E2.a4/E2.a3; E3 = ellchangecurve(E2,[1,0,s,0]);
  u = E3.a3/E3.a2; E4 = ellchangecurve(E3,[u,0,0,0]);
  [-E4.a2, 1-E4.a1]
};
\\ all points of exact order N (N = torsN(f)) of E(Q)
toppts(E,N) = {
  my(T = elltors(E), G = T[2], gens = T[3], pts = List());
  if(#G==1, for(k=1,G[1]-1, my(R=ellmul(E,gens[1],k)); if(R!=[0] && ellorder(E,R)==N, listput(pts,R))),
    for(i=0,G[1]-1, for(j=0,G[2]-1, my(R=elladd(E, ellmul(E,gens[1],i), ellmul(E,gens[2],j))); if(R!=[0] && ellorder(E,R)==N, listput(pts,R)))));
  Vec(pts)
};
\\ all family parameters x (as rationals, oo allowed as [1,0]) whose fibre is isomorphic to E
params(f,E) = {
  my(N=torsN(f), res=List(), BC=tbc(f,'x));
  if(elltors(E)[2] != torsG(f), return([]));
  foreach(toppts(E,N), R, my(bc = tate(E,R), g);
    g = gcd(numerator(BC[2]-bc[2]), numerator(BC[1]-bc[1]));
    foreach(nfroots(, g), r, my(bb=tbc(f,r)); if(bb==bc, listput(res, r))));
  Set(Vec(res))
};
\\ ---- symmetries: Moebius maps x -> (al*x+be)/(ga*x+de) fixing the fibre (found by symfind.gp, checked) ----
SYMS = [[[0, 1, -1, 1], [1, -1, 1, 0], [1, 0, 0, 1]], [[-1, 1, 0, 1], [1, 0, 0, 1]], [[0, 1, -1, 1], [1, -1, 1, 0], [1, 0, 0, 1]], [[-1, 1, -2, 1], [1, 0, 0, 1]], [[-1, 1, 0, 1], [1, 0, 0, 1]], [[-5, 9, -1, 5], [-5, 21, -1, 1], [-1, -15, -1, 1], [-1, 6, 0, 1], [-1, 21, -1, 5], [1, 0, 0, 1]], [[-4, -1, 8, 0], [-4, -1, 8, 4], [-2, -1, 0, 2], [-2, -1, 8, 2], [-1, 0, 4, 1], [0, -1, 8, 4], [0, 1, 8, 0], [1, 0, 0, 1]]];  \\ from symfind.gp; re-checked by selftest.gp (params(f,E) == symorbit)
syms(f) = SYMS[famidx(f)];
mob(M,P,Q) = { my(a=M[1]*P+M[2]*Q, b=M[3]*P+M[4]*Q, g=gcd(a,b)); if(g==0, return([0,0])); a/=g; b/=g; if(b<0 || (b==0 && a<0), a=-a; b=-b); [a,b] };
symorbit(f,P,Q) = Set(vector(#syms(f), i, mob(syms(f)[i],P,Q)));
\\ canonical representative: smallest max(|p|,|q|), ties broken by [q,p] lexicographic
canon(f,P,Q) = { my(o=symorbit(f,P,Q), best=o[1]); foreach(o, v, my(m=max(abs(v[1]),abs(v[2])), mb=max(abs(best[1]),abs(best[2])));
  if(m<mb || (m==mb && lex([v[2],v[1]], [best[2],best[1]]) < 0), best=v)); best };
