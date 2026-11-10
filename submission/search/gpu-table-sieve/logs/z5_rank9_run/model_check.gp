\\ Integral model of the sieve family at the found parameters, minimal model, and discriminant factorisation built from the factors p, q, p+q, p+9q (resp. p^2-11pq-q^2) (usage: FAM=66|5 P=<p> Q=<q> ECFILE=<file> gp -q model_check.gp)
default(parisizemax, 2*10^9);
fam = eval(getenv("FAM")); p = eval(getenv("P")); q = eval(getenv("Q"));
L = readstr(getenv("ECFILE")); sub = 0;
for(i=1,#L, my(v=strsplit(L[i]," ")); if(#L[i]>1 && v[1]=="E", sub=[eval(v[2]),eval(v[3]),eval(v[4]),eval(v[5]),eval(v[6])]));
if(fam==66, a=[p-q, -q*(p+q), -p*q*(p+q), 0, 0]; extra = p+9*q, a=[q-p, -p*q, -p*q^2, 0, 0]; extra = p^2-11*p*q-q^2);
E = ellinit(a);
print("integral model [a1,a2,a3,a4,a6]=", a);
Em = ellminimalmodel(E, &w);
print("minimal model ", Em[1..5]);
print("same as submitted: ", Em[1..5]==sub);
print("change of variables [u,r,s,t] integral -> minimal: ", w);
Hn(E)=log(max(abs(E.c4)^3, E.c6^2))/12.;
print("Hint=", Hn(E));
print("Hmin=", Hn(Em));
cands = Set(concat([[2,3,5,7], factor(abs(p))[,1]~, factor(abs(q))[,1]~, factor(abs(p+q))[,1]~, factor(abs(extra))[,1]~]));
fac(D)={ my(res=List(), r=abs(D), v); for(i=1,#cands, v=valuation(r,cands[i]); if(v>0, listput(res,[cands[i],v]); r/=cands[i]^v)); [Mat(Vec(res)~), r] };
\\ Mat of a list of pairs
fmat(D)={ my(res=List(), r=abs(D), v, M); for(i=1,#cands, v=valuation(r,cands[i]); if(v>0, listput(res,[cands[i],v]); r/=cands[i]^v)); M=matrix(#res,2,i,j,res[i][j]); [M, r] };
F0 = fmat(E.disc); Fm = fmat(Em.disc);
print("sign of the discriminant: ", sign(E.disc), "   (the factorisations below are of the absolute value)");
print("Delta0 = ", F0[1]~, "   (cofactor left after dividing out all primes of p, q, p+q, the last factor and 2,3,5,7: ", F0[2], ")");
print("Delta_min = ", Fm[1]~, "   (cofactor: ", Fm[2], ")");
print("p = ", factor(abs(p)), "   q = ", factor(abs(q)));
print("p+q = ", sign(p+q), " * ", factor(abs(p+q)), "   ", if(fam==66,"p+9q = ","p^2-11pq-q^2 = "), sign(extra), " * ", factor(abs(extra)));
T = elltors(Em); print("torsion order: ", T[1], " structure ", T[2]);
bad = 0; forprime(l=2,1000, if(Em.disc%l==0, bad++));
print("good primes among the first 168 (<1000): ", 168-bad);
