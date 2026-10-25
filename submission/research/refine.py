import ast,collections,json,pathlib,sys
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[2]))
from ec_verify import gp,records
ROOT=pathlib.Path("submission/research")
buckets=collections.defaultdict(list)
table=records.load()
for line in (ROOT/"sieve-wide.txt").read_text().splitlines():
 if not line.startswith("["):continue
 row=ast.literal_eval(line)
 f,p,q,h,s,a=row
 target=next(r for r in range(1,10) if records.baseline(table,f"{f}/{r}") is None or round(h,6)<records.baseline(table,f"{f}/{r}"))
 buckets[f,target].append(row)
selected=[]
for key,group in buckets.items():
 selected+=sorted(group,key=lambda x:x[4])[:500]
body=['ps=primes([2,5000]);','ws=vector(#ps,i,log(ps[i])/ps[i]);']
for i,row in enumerate(selected):
 body+= [f'E=ellinit({row[5]});s=0;for(i=1,#ps,if(E.disc%ps[i],s+=ellap(E,ps[i])*ws[i]));print("S{i}=",s);']
out=gp.run_gp("\n".join(body),timeout=180)
with (ROOT/"sieve-refined.txt").open("w") as stream:
 for i,row in enumerate(selected):
  row[4]=float(out[f"S{i}"])
  stream.write(str(row)+"\n")
print("REFINED",len(selected),flush=True)
