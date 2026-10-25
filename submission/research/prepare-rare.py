import ast,collections,pathlib,sys
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[2]))
from ec_verify import gp,records
root=pathlib.Path("submission/research")
buckets=collections.defaultdict(list)
for f in (7,10,12):
 for line in (root/f"cached-rare-{f}.txt").read_text().splitlines():
  f,p,q,h,s=line.split();f,p,q=int(f),int(p),int(q);h,s=float(h),float(s)
  target=next(r for r in range(1,8) if records.baseline(records.load(),f"{f}/{r}") is None or round(h,6)<records.baseline(records.load(),f"{f}/{r}"))
  buckets[f,target].append((f,p,q,h,s))
rows=[]
for key,group in buckets.items():
 rows+=sorted(group,key=lambda x:x[4])[:600]
print("REFINING_RARE",len(rows),flush=True)
body=['ps=primes([2,5000]);','ws=vector(#ps,i,log(ps[i])/ps[i]);']
def model(f,p,q):
 if f==7:
  a2=-p*p*(p-q)*q
  return [q*q-p*p+p*q,a2,a2*q*q,0,0]
 if f==10:
  d=p*q-(p-q)**2;nc=p*(p-q)*(2*p-q);a2=-nc*p*p*q
  return [q*d-nc,a2,a2*q*d,0,0]
 d=p-q;nc=p*(q-2*p)*(3*p*p-3*p*q+q*q);nb=nc*(2*p*q-2*p*p-q*q);a2=-nb*d*d
 return [q*d**3-nc,a2,a2*q*d**3,0,0]
for i,(f,p,q,h,s) in enumerate(rows):
 a=model(f,p,q)
 body += [f'E=ellminimalmodel(ellinit({a}));h=log(max(abs(E.c4)^3,E.c6^2))/12;s=0;for(i=1,#ps,if(E.disc%ps[i],s+=ellap(E,ps[i])*ws[i]));print("C{i}=",[h,s,[E.a1,E.a2,E.a3,E.a4,E.a6]]);']
out=gp.run_gp("\n".join(body),timeout=240)
with (root/"sieve-cached-rare.txt").open("w") as stream:
 for i,(f,p,q,h,s) in enumerate(rows):
  hm,score,a=ast.literal_eval(out[f"C{i}"])
  stream.write(str([f,p,q,hm,score,a])+"\n")
print("REFINED_RARE",len(rows),flush=True)
checks=rows[::max(1,len(rows)//20)]
body=['ps=primes([2,1000]);','ws=vector(#ps,i,log(ps[i])/ps[i]);']
for i,(f,p,q,h,s) in enumerate(checks):
 body += [f'E=ellinit({model(f,p,q)});s=0;for(i=1,#ps,if(E.disc%ps[i],s+=ellap(E,ps[i])*ws[i]));print("S{i}=",s);']
out=gp.run_gp("\n".join(body),timeout=20)
error=max(abs(float(out[f"S{i}"])-row[4]) for i,row in enumerate(checks))
print("RARE_CACHE_DIRECT_GP_MAX_ERROR",error,flush=True)
assert error<1e-10
