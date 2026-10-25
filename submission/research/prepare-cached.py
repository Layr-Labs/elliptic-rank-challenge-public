"""Normalize the C++ shortlist and independently recompute its prime sums in GP."""
import ast,pathlib,sys
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[2]))
from ec_verify import gp
root=pathlib.Path("submission/research")
rows=[]
for line in (root/"cached-shortlist.txt").read_text().splitlines():
 f,p,q,h,s=line.split();rows.append((int(f),int(p),int(q),float(h),float(s)))
body=['ps=primes([2,5000]);','ws=vector(#ps,i,log(ps[i])/ps[i]);']
for i,(f,p,q,h,s) in enumerate(rows):
 a=[q-p,-p*q,-p*q*q,0,0] if f==5 else [q-p,-p*(p+q),-p*q*(p+q),0,0]
 body+= [f'E=ellminimalmodel(ellinit({a}));h=log(max(abs(E.c4)^3,E.c6^2))/12;s=0;for(i=1,#ps,if(E.disc%ps[i],s+=ellap(E,ps[i])*ws[i]));print("C{i}=",[h,s,[E.a1,E.a2,E.a3,E.a4,E.a6]]);']
out=gp.run_gp("\n".join(body),timeout=240)
with (root/"sieve-cached.txt").open("w") as stream:
 for i,(f,p,q,h,s) in enumerate(rows):
  hm,score,a=ast.literal_eval(out[f"C{i}"])
  stream.write(str([f,p,q,hm,score,a])+"\n")
print("NORMALIZED_AND_REFINED",len(rows),flush=True)
# Compare cached S(1000) against direct GP arithmetic on a reproducible sample.
checks=rows[::max(1,len(rows)//20)]
body=['ps=primes([2,1000]);','ws=vector(#ps,i,log(ps[i])/ps[i]);']
for i,(f,p,q,h,s) in enumerate(checks):
 a=[q-p,-p*q,-p*q*q,0,0] if f==5 else [q-p,-p*(p+q),-p*q*(p+q),0,0]
 body += [f'E=ellinit({a});s=0;for(i=1,#ps,if(E.disc%ps[i],s+=ellap(E,ps[i])*ws[i]));print("S{i}=",s);']
out=gp.run_gp("\n".join(body),timeout=20)
error=max(abs(float(out[f"S{i}"])-row[4]) for i,row in enumerate(checks))
print("CACHE_DIRECT_GP_MAX_ERROR",error,flush=True)
assert error<1e-10
