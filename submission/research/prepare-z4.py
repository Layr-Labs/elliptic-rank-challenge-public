import ast,collections,json,pathlib,sys
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[2]))
from ec_verify import gp,records
root=pathlib.Path("submission/research")
parent=json.loads((root/"initial-board.json").read_text())
thresholds={x["cell"]:x["naive_height"] for x in parent["metrics"]["curves"] if x["contribution"]>0}
groups=collections.defaultdict(list)
for line in (root/"cached-z4.txt").read_text().splitlines():
 i,p,q,h,s=line.split();i,p,q=int(i),int(p),int(q);h,s=float(h),float(s)
 target=next(r for r in range(1,20) if records.baseline(records.load(),f"4/{r}") is None or h<min(records.baseline(records.load(),f"4/{r}"),thresholds.get(f"4/{r}",float("inf"))))
 groups[i,target].append((i,p,q,h,s))
rows=[]
for key,group in groups.items():
 rows+=sorted(group,key=lambda x:x[4])[:600]
print("REFINING_Z4",len(rows),flush=True)
body=['ps=primes([2,5000]);ws=vector(#ps,i,log(ps[i])/ps[i]);']
for j,(i,p,q,h,s) in enumerate(rows):
 if i==1:
  aa=(8*p-q)*(32*p+7*q);bb=8*q*(p+q)*(15*p-8*q)*(31*p-7*q)
 else:
  aa=-8*q*(80*p+9*q);bb=2*(p-2*q)*(2*p-q)*(18*p-q)*(2*p-81*q)
 a=[aa,bb,aa*bb,0,0]
 body += [f'E=ellminimalmodel(ellinit({a}));h=log(max(abs(E.c4)^3,E.c6^2))/12;s=0;for(ii=1,#ps,if(E.disc%ps[ii],s+=ellap(E,ps[ii])*ws[ii],if(ps[ii]<=1000,s-=0.4)));print("C{j}=",[h,s,[E.a1,E.a2,E.a3,E.a4,E.a6]]);']
out=gp.run_gp("\n".join(body),timeout=300)
seen=set()
with (root/"sieve-z4.txt").open("w") as stream:
 for j,(i,p,q,h,s) in enumerate(rows):
  hm,score,a=ast.literal_eval(out[f"C{j}"])
  if tuple(a) in seen:continue
  seen.add(tuple(a))
  stream.write(str([4,p,q,hm,score,a,f"z4:{i}",0])+"\n")
print("REFINED_Z4",len(seen),flush=True)
