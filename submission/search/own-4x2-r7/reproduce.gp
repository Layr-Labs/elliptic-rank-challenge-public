\\ Reproduce our independently searched fiber; stdout is a candidate .ec.
\\ Published family: y^2=x(x+p^2)(x+q^2).
\\ Run from repository root with pinned PARI/GP 2.17.4.
p=296037;q=1024556;
E=ellinit([0,p^2+q^2,0,p^2*q^2,0]);
M=ellminimalmodel(E,&ch);
if(M[1..5]!=[1,0,0,-21199732376797362775114,1157881752110224304095839182626259],error("minimal model mismatch"));
if(elltors(M)[2]!=[4,2],error("torsion mismatch"));
setrand(1);R=ellrank(M,1);
if(R[1]!=7||R[2]!=7||#R[4]!=7,error("rank search did not reproduce seven points"));
if(!vecprod(vector(7,k,ellisoncurve(M,R[4][k]))),error("point off curve"));
print("# Independently searched parameter: family=44 p=296037 q=1024556");
print("E ",strjoin(apply(x->Str(x),M[1..5])," "));
for(k=1,7,print("P ",R[4][k][1]," ",R[4][k][2]));
print("#CLAIM rank=7");
quit;
