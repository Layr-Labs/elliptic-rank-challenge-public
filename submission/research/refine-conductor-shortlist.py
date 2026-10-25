"""Refine both SIMD ranking modes; never treat a heuristic as a certificate."""
import ast,collections,concurrent.futures,json,pathlib,sys,time,hashlib
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[2]))
from ec_verify import gp,records
root=pathlib.Path("submission/research")
parent=json.loads((root/"initial-board.json").read_text())
thresholds={r["cell"]:r["naive_height"] for r in parent["metrics"]["curves"] if r["contribution"]>0}
unique={}
for line in (root/"full4-conductor-shortlist.txt").read_text().splitlines():
 _,p,q,h,s,lnN,mode=line.split();p,q=int(p),int(q)
 unique.setdefault((p,q),(p,q,float(h),float(s),float(lnN)))
rows=list(unique.values())
print("UNIQUE_CONDUCTOR_SHORTLIST",len(rows),flush=True)
refine_cache=pathlib.Path("/tmp/ec-rank-conductor-refine-cache")
refine_cache.mkdir(exist_ok=True)
def refine(batch):
 cache_key=hashlib.sha256(json.dumps(batch,separators=(",",":")).encode()).hexdigest()
 cache_file=refine_cache/(cache_key+".json")
 if cache_file.exists():return json.loads(cache_file.read_text())
 body=['ps=primes([2,10000]);ws=vector(#ps,i,log(ps[i])/ps[i]);']
 for j,(p,q,h,s,lnN) in enumerate(batch):
  a=[0,2*(p*p+q*q),0,(p*p-q*q)**2,0]
  body += [f'E=ellminimalmodel(ellinit({a}));h=log(max(abs(E.c4)^3,E.c6^2))/12;s=0;nb=0;for(ii=1,#ps,if(E.disc%ps[ii],s+=ellap(E,ps[ii])*ws[ii],if(ps[ii]<=2000,nb++)));print("C{j}=",[h,s,nb,[E.a1,E.a2,E.a3,E.a4,E.a6]]);']
 out=gp.run_gp("\n".join(body),timeout=120)
 result=[(raw,ast.literal_eval(out[f"C{j}"])) for j,raw in enumerate(batch)]
 cache_file.write_text(json.dumps(result,separators=(",",":")))
 return result
allrows=[];batches=[rows[i:i+1000] for i in range(0,len(rows),1000)]
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
 futures=[pool.submit(refine,batch) for batch in batches]
 for future in concurrent.futures.as_completed(futures):
  allrows+=future.result()
  print("REFINEMENT_PROGRESS",len(allrows),"/",len(rows),flush=True)
groups=collections.defaultdict(list);by_model={}
for raw,result in allrows:
 p,q,_,s2000,lnN=raw;h,s10000,nb,a=result
 target=next(r for r in range(1,10) if records.baseline(records.load(),f"4x2/{r}") is None or round(h,6)<min(records.baseline(records.load(),f"4x2/{r}"),thresholds.get(f"4x2/{r}",float("inf"))))
 record=dict(p=p,q=q,h=h,s10000=s10000,small_bad=nb,log_conductor=lnN,ainvs=a,target=target)
 key=tuple(a)
 # Exact model deduplication only; full board fingerprints handle isogenies/twists.
 if key not in by_model:by_model[key]=record
for record in by_model.values():groups[record["target"]].append(record)
selected={}
for target,group in groups.items():
 for mode in (0,1):
  score=lambda row:row["s10000"]-(0.4*row["small_bad"] if mode==0 else 0.5*row["log_conductor"])
  for row in sorted(group,key=score)[:3000]:
   selected[tuple(row["ainvs"])]=row
with (root/"sieve-full4-conductor.txt").open("w") as stream:
 for row in selected.values():
  priority=row["s10000"]-0.5*row["log_conductor"]
  stream.write(str(["4x2",row["p"],row["q"],row["h"],priority,row["ainvs"]])+"\n")
(root/"conductor-refined-metrics.json").write_text(json.dumps(list(selected.values()),separators=(",",":")))
print("SELECTED_UNION",len(selected),"ALL_UNIQUE_MODELS",len(by_model),flush=True)
