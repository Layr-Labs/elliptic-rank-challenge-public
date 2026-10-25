import ast,collections,pathlib,sys
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[2]))
from ec_verify import gp,records
root=pathlib.Path("submission/research")
buckets=collections.defaultdict(list)
for line in (root/"cached-full4.txt").read_text().splitlines():
 _,p,q,h,s=line.split();p,q=int(p),int(q);h,s=float(h),float(s)
 target=next(r for r in range(1,10) if records.baseline(records.load(),f"4x2/{r}") is None or round(h,6)<records.baseline(records.load(),f"4x2/{r}"))
 buckets[target].append((p,q,h,s))
rows=[]
for key,group in buckets.items():
 rows+=sorted(group,key=lambda x:x[3])[:1000]
print("REFINING_FULL4",len(rows),flush=True)
body=['ps=primes([2,5000]);','ws=vector(#ps,i,log(ps[i])/ps[i]);']
def model(p,q):
 return [0,2*(p*p+q*q),0,(p*p-q*q)**2,0]
for i,(p,q,h,s) in enumerate(rows):
 body += [f'E=ellminimalmodel(ellinit({model(p,q)}));h=log(max(abs(E.c4)^3,E.c6^2))/12;s=0;for(i=1,#ps,if(E.disc%ps[i],s+=ellap(E,ps[i])*ws[i]));print("C{i}=",[h,s,[E.a1,E.a2,E.a3,E.a4,E.a6]]);']
out=gp.run_gp("\n".join(body),timeout=240)
with (root/"sieve-full4.txt").open("w") as stream:
 seen=set()
 for i,(p,q,h,s) in enumerate(rows):
  hm,score,a=ast.literal_eval(out[f"C{i}"])
  if tuple(a) in seen:continue
  seen.add(tuple(a))
  stream.write(str(["4x2",p,q,hm,score,a])+"\n")
print("REFINED_FULL4",len(seen),flush=True)
checks=rows[::max(1,len(rows)//20)]
body=['ps=primes([2,1000]);','ws=vector(#ps,i,log(ps[i])/ps[i]);']
for i,(p,q,h,s) in enumerate(checks):
 body += [f'E=ellinit({model(p,q)});s=0;for(i=1,#ps,if(E.disc%ps[i],s+=ellap(E,ps[i])*ws[i]));print("S{i}=",s);']
out=gp.run_gp("\n".join(body),timeout=20)
error=max(abs(float(out[f"S{i}"])-row[3]) for i,row in enumerate(checks))
print("FULL4_CACHE_DIRECT_GP_MAX_ERROR",error,flush=True)
assert error<1e-10
