\\ famsubst.gp - generate the fams.txt blocks of a family re-expressed in coordinates (pk, qk) with qk = alpha*p + beta*q (an "outer factor" variant), pk completing it to a unimodular basis.
\\ usage: gp -q famsubst.gp   with  FAMBASE, NEWID, ALPHA, BETA set on the command line via -D? (here: edit the call list at the bottom)
fams = readstr("fams.txt");
getfam(F) = { my(i=1, A=vector(3), FD=0); while(i<=#fams, my(t=strsplit(fams[i]," ")); if(t[1]=="FAM" && eval(t[2])==F, FD=eval(t[3]); for(k=1,3, my(s=select(x->#x>0,strsplit(fams[i+k]," ")), n=eval(s[2])); A[k]=vector(n,j,eval(s[2+j]))); break); i++); [FD,A] }
\\ new family: returns the text of the FAM block
variant(F, newid, al, be) = {
  my(r=getfam(F), FD=r[1], A=r[2], g, ga, de, w, P, Q, pk, qk, out);
  \\ find (ga,de) with al*de - be*ga = 1
  [ga, de] = [0,0]; my(bz = bezout(al, be)); \\ bezout gives [u,v,g] with u*al+v*be=g
  my(u=bz[1], v=bz[2], gg=bz[3]); if(gg != 1 && gg != -1, error("alpha,beta not coprime")); u*=gg; v*=gg;   \\ u*al + v*be = 1
  \\ new variables: qk = al*p + be*q ; pk = -v*p + u*q  (det [[-v,u],[al,be]] = -v*be - u*al = -1)
  \\ inverse: p = -be*pk... solve with matrix
  my(M = [-v, u; al, be], Mi = M^(-1));        \\ (pk,qk)^T = M (p,q)^T  ->  (p,q)^T = Mi (pk,qk)^T
  pk = 'pk; qk = 'qk;
  my(pp = Mi[1,1]*pk + Mi[1,2]*qk, qq = Mi[2,1]*pk + Mi[2,2]*qk);
  out = Str("FAM ", newid, " ", FD, "\n");
  for(k=1,3, w = k*FD;
    my(poly = sum(j=1,#A[k], A[k][j]*pp^(j-1)*qq^(w-(j-1))));
    my(cv = vector(w+1, j, polcoef(polcoef(poly, j-1, 'pk), w-(j-1), 'qk)));
    out = Str(out, "A", k, " ", w+1, " ", strjoin(apply(x->Str(x), cv), " "), "\n"));
  out
}
\\ self-test against the hand-made FAM 66 (p) and FAM 67 (p+q) of Z/6:
print(variant(6, 66, 1, 0));
print(variant(6, 67, 1, 1));
print(variant(26, 261, 1, -1));
print(variant(26, 265, 1, -5));
