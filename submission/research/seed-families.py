import ast,json,pathlib,sys
from fractions import Fraction
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[2]))
from ec_verify import gp
root=pathlib.Path("submission/research")
rows=[]
for row in map(ast.literal_eval,(root/"z3-families.txt").open()):
 if row[-1]<3:continue
 den=Fraction(row[3]).denominator
 m=next((m for m in range(den,301,den) if (m**3*Fraction(row[4])).denominator==1),None)
 if m is not None:rows.append((row,m))
body=['read("/home/linux/ec-rank/submission/research/z3-family.gp");']
for i,(row,m) in enumerate(rows):
 r,s,v,lam,k,count=row
 body += [f'F=z3family({r},{s},{v},{lam},{k},101/103);E=ellminimalmodel(F[1],&vc);P=vector(#F[2],ii,ellchangepoint(F[2][ii],vc));ev=qfjacobi(ellheightmatrix(E,P))[1];rankseed=sum(ii=1,#ev,ev[ii]>10^-20);if(rankseed>=4,F=z3family({r},{s},{v},{lam},{k},101);E=ellminimalmodel(F[1]);h=log(max(abs(E.c4)^3,E.c6^2))/12;print("F{i}=",[rankseed,h-2*log(101)]));']
out=gp.run_gp("\n".join(body),timeout=180,precision=57)
selected=[]
for i,(row,m) in enumerate(rows):
 if f"F{i}" in out:
  rank,offset=ast.literal_eval(out[f"F{i}"])
  selected.append(dict(roots=row[:3],lam=row[3],k=row[4],seed_rank=rank,m=m,offset=offset))
selected.sort(key=lambda x:(-x["seed_rank"],x["offset"]))
(root/"seed-ranked-families.json").write_text(json.dumps(selected,indent=1))
print("SEED_RANKS",__import__("collections").Counter(x["seed_rank"] for x in selected),flush=True)
print(json.dumps(selected[:15],indent=1),flush=True)
