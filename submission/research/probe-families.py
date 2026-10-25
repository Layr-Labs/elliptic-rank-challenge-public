import ast,json,pathlib,sys
from fractions import Fraction
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[2]))
from ec_verify import gp
root=pathlib.Path("submission/research")
rows=[r for r in map(ast.literal_eval,(root/"z3-families.txt").open()) if r[-1]>=3]
if "--small" in sys.argv:
 rows=[row for row in rows if any((m**3*Fraction(row[4])).denominator==1 for m in range(Fraction(row[3]).denominator,61,Fraction(row[3]).denominator))]
body=['read("/home/linux/ec-rank/submission/research/z3-family.gp");']
for i,(r,s,v,lam,k,count) in enumerate(rows):
 body += [f"""lower=99;
for(jj=1,2,tprobe=if(jj==1,2/3,3/7);F=z3family({r},{s},{v},{lam},{k},tprobe);if(#F[1]==0,lower=0;break);E=ellminimalmodel(F[1],&vc);P=vector(#F[2],ii,ellchangepoint(F[2][ii],vc));R=ellrank(E,0,P);lower=min(lower,R[1]));
if(lower>=5,m=lcm(denominator({lam}),prod(jj=1,matsize(factor(denominator({k})))[1],factor(denominator({k}))[jj,1]^ceil(factor(denominator({k}))[jj,2]/3)));F=z3family({r},{s},{v},{lam},{k},101);E=ellminimalmodel(F[1]);h=log(max(abs(E.c4)^3,E.c6^2))/12;print("F{i}=",[lower,m,h-2*log(101)]));
"""]
out=gp.run_gp("\n".join(body),timeout=240)
selected=[]
for i,row in enumerate(rows):
 if f"F{i}" in out:
  rank,m,offset=ast.literal_eval(out[f"F{i}"])
  selected.append(dict(roots=row[:3],lam=row[3],k=row[4],rank=rank,m=m,offset=offset))
selected.sort(key=lambda x:x["offset"])
(root/("constructed-families-small.json" if "--small" in sys.argv else "constructed-families.json")).write_text(json.dumps(selected,indent=1))
print("FAMILIES_WITH_TWO_RANK5_SPECIALIZATIONS",len(selected),flush=True)
print(json.dumps(selected[:12],indent=1),flush=True)
