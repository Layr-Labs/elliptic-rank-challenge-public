import ast,collections,json,pathlib,sys
from fractions import Fraction
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[2]))
from ec_verify import gp,records
root=pathlib.Path("submission/research")
families=json.loads((root/"constructed-families-small.json").read_text())
groups=collections.defaultdict(list)
for line in (root/"cached-constructed.txt").read_text().splitlines():
 i,p,q,h,s=line.split();i,p,q=int(i),int(p),int(q);h,s=float(h),float(s)
 target=next(r for r in range(1,20) if records.baseline(records.load(),f"3/{r}") is None or round(h,6)<records.baseline(records.load(),f"3/{r}"))
 groups[i,target].append((i,p,q,h,s))
rows=[]
for key,group in groups.items():
 rows+=sorted(group,key=lambda x:x[4])[:350]
print("REFINING_CONSTRUCTED",len(rows),flush=True)
body=['ps=primes([2,5000]);ws=vector(#ps,i,log(ps[i])/ps[i]);']
for j,(i,p,q,h,s) in enumerate(rows):
 f=families[i-1];r1,r2,r3=f["roots"];m=f["m"];b=int(m*Fraction(f["lam"]));k=int(m**3*Fraction(f["k"]))
 a=[m*p*p+b*q*q,0,k*(p*p-r1*r1*q*q)*(p*p-r2*r2*q*q)*(p*p-r3*r3*q*q),0,0]
 body += [f'E=ellminimalmodel(ellinit({a}));h=log(max(abs(E.c4)^3,E.c6^2))/12;s=0;for(ii=1,#ps,if(E.disc%ps[ii],s+=ellap(E,ps[ii])*ws[ii]));print("C{j}=",[h,s,[E.a1,E.a2,E.a3,E.a4,E.a6]]);']
out=gp.run_gp("\n".join(body),timeout=300)
seen=set()
with (root/"sieve-constructed-wide.txt").open("w") as stream:
 for j,(i,p,q,h,s) in enumerate(rows):
  hm,score,a=ast.literal_eval(out[f"C{j}"])
  if tuple(a) in seen:continue
  seen.add(tuple(a));f=families[i-1];r1,r2,r3=f["roots"]
  meta=f'z3:{r1},{r2},{r3},{f["lam"]},{f["k"]}'
  stream.write(str([3,p,q,hm,score,a,meta,0])+"\n")
print("REFINED_CONSTRUCTED",len(seen),flush=True)
