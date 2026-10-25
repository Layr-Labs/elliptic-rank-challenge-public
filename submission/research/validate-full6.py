import json,math,pathlib,random,sys
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[2]))
from ec_verify import gp
lookup={}
for line in open("/tmp/ec-rank-full6-local-tables.txt"):
 r,s,u,f=map(int,line.split());lookup[r,s]=(u,f)
rng=random.Random(20261025);samples=[]
while len(samples)<100:
 p=rng.randint(-50000,50000);q=rng.randint(1,50000)
 if math.gcd(p,q)!=1 or p==0 or any(x==0 for x in [p-q,p+q,3*p-5*q,p-3*q,3*p-q]):continue
 samples.append((p,q))
body=['read("/home/linux/ec-rank/submission/research/full6-family.gp");']
for i,(p,q) in enumerate(samples):
 body.append(f'E=full6curve({p},{q});red2=elllocalred(E,2);print("U_{i}=",valuation(red2[3][1],2));print("F2_{i}=",red2[1]);print("F3_{i}=",elllocalred(E,3)[1]);Em=ellminimalmodel(E);print("H_{i}=",log(max(abs(Em.c4)^3,Em.c6^2))/12);print("LOGN_{i}=",log(ellglobalred(E)[1]));')
out=gp.run_gp("\n".join(body),timeout=60)
def radlog(n):
 n=abs(n);total=0.
 for prime in (2,3):
  while n%prime==0:n//=prime
 d=5
 while d*d<=n:
  if n%d==0:
   total+=math.log(d)
   while n%d==0:n//=d
  d+=2
 if n>1:total+=math.log(n)
 return total
height_err=[];conductor_err=[]
for i,(p,q) in enumerate(samples):
 A=37*q**4-84*p*q**3+102*p*p*q*q-36*p**3*q-3*p**4
 B=32*q*(p-q)**3*(p+q)**3*(3*p-5*q)
 c4=16*(A*A-3*B);c6=-64*A**3+288*A*B
 u,f=lookup[p%128,q%128];f3=2 if q%3==0 else 1
 h=max(math.log(abs(c4))/4,math.log(abs(c6))/6)-u*math.log(2)
 height_err.append(abs(h-float(out[f"H_{i}"])))
 assert u==int(out[f"U_{i}"]) and f==int(out[f"F2_{i}"]) and f3==int(out[f"F3_{i}"])
 lnN=sum(radlog(n) for n in [q,p-q,p+q,3*p-5*q,p-3*q,3*p-q])+f*math.log(2)+f3*math.log(3)
 conductor_err.append(abs(lnN-float(out[f"LOGN_{i}"])))
assert max(height_err)<1e-10 and max(conductor_err)<1e-10
report=dict(samples=100,max_height_error=max(height_err),max_log_conductor_error=max(conductor_err))
pathlib.Path("submission/research/full6-validation.json").write_text(json.dumps(report,indent=1))
print("FULL6_VALIDATED",report,flush=True)
