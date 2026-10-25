import ast,json,math,pathlib,sys
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[2]))
from ec_verify import gp
tables=[]
for line in open("/tmp/ec-rank-full6-prime-tables.txt"):
 p,v=line.split(" ",1);tables.append((int(p),ast.literal_eval(v)))
rows=[s.split() for s in open("/tmp/ec-rank-full6-validation-shortlist.txt")]
selected={}
for sign in (-1,1):
 matching=[row for row in rows if int(row[1])*sign>0]
 for row in matching[::max(1,len(matching)//25)][:25]:selected[int(row[1]),int(row[2])]=row
assert any(p<0 for p,q in selected) and any(p>0 for p,q in selected)
max_counter_error=0
for (p,q),row in selected.items():
 total=0
 for prime,traces in tables:
  if q%prime==0:continue
  trace=traces[(p*pow(q,-1,prime))%prime]
  if trace==9999:continue
  value=128*trace*math.log(prime)/prime
  total+=math.floor(value+0.5) if value>=0 else math.ceil(value-0.5)
 actual=round(float(row[4])*128)
 max_counter_error=max(max_counter_error,abs(total-actual))
assert max_counter_error==0
sample=list(selected)[:10]+[pair for pair in selected if pair[0]>0][:5]
body=['read("/home/linux/ec-rank/submission/research/full6-family.gp");']
for i,(p,q) in enumerate(sample):
 body.append(f'E=ellminimalmodel(full6curve({p},{q}));s=0;forprime(pr=2,2000,if(E.disc%pr,s+=ellap(E,pr)*log(pr)/pr));print("S_{i}=",s);')
out=gp.run_gp("\n".join(body),timeout=30)
errors=[]
for i,(p,q) in enumerate(sample):
 value=0.
 for prime,traces in tables:
  if q%prime==0:continue
  trace=traces[(p*pow(q,-1,prime))%prime]
  if trace!=9999:value+=trace*math.log(prime)/prime
 errors.append(abs(value-float(out[f"S_{i}"])))
assert max(errors)<1e-10
report=dict(counter_samples=len(selected),negative_samples=sum(p<0 for p,q in selected),max_integer_counter_error=max_counter_error,gp_samples=len(sample),max_gp_trace_error=max(errors))
pathlib.Path("submission/research/full6-simd-validation.json").write_text(json.dumps(report,indent=1))
print("FULL6_SIMD_VALIDATED",report,flush=True)
