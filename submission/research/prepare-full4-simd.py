import ast,collections,json,pathlib,sys
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[2]))
from ec_verify import gp,records
root=pathlib.Path("submission/research")
parent=json.loads((root/"initial-board.json").read_text())
thresholds={r["cell"]:r["naive_height"] for r in parent["metrics"]["curves"] if r["contribution"]>0}
groups=collections.defaultdict(list)
for line in (root/"full4-simd-shortlist.txt").read_text().splitlines():
 _,p,q,h,s=line.split();p,q=int(p),int(q);h,s=float(h),float(s)
 target=next(r for r in range(1,10) if records.baseline(records.load(),f"4x2/{r}") is None or round(h,6)<min(records.baseline(records.load(),f"4x2/{r}"),thresholds.get(f"4x2/{r}",float("inf"))))
 groups[target].append((p,q,h,s))
rows=[]
for target,group in groups.items():rows+=sorted(group,key=lambda x:x[3])[:4000]
print("REFINING_SIMD",len(rows),flush=True)
body=['ps=primes([2,5000]);ws=vector(#ps,i,log(ps[i])/ps[i]);']
for j,(p,q,h,s) in enumerate(rows):
 a=[0,2*(p*p+q*q),0,(p*p-q*q)**2,0]
 body += [f'E=ellminimalmodel(ellinit({a}));h=log(max(abs(E.c4)^3,E.c6^2))/12;s=0;for(ii=1,#ps,if(E.disc%ps[ii],s+=ellap(E,ps[ii])*ws[ii],if(ps[ii]<=2000,s-=0.4)));print("C{j}=",[h,s,[E.a1,E.a2,E.a3,E.a4,E.a6]]);']
out=gp.run_gp("\n".join(body),timeout=300)
seen=set()
with (root/"sieve-full4-simd.txt").open("w") as stream:
 for j,(p,q,h,s) in enumerate(rows):
  hm,score,a=ast.literal_eval(out[f"C{j}"])
  if tuple(a) in seen:continue
  seen.add(tuple(a));stream.write(str(["4x2",p,q,hm,score,a])+"\n")
print("REFINED_SIMD",len(seen),flush=True)
