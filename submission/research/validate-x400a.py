import ast,json,math,pathlib,random,sys
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[2]))
from ec_verify import gp
rows=[s.split() for s in pathlib.Path('submission/research/x400a-shortlist.txt').read_text().splitlines()]
random.seed(20261025)
selected=random.sample(rows,30)
tables={}
for line in pathlib.Path('/tmp/ec-rank-x400a-prime-tables.txt').read_text().splitlines():
 pr,v=line.split(' ',1);pr=int(pr);tables[pr]=ast.literal_eval(v)
body=['ps=primes([2,2000]);']
checks=[]
for i,row in enumerate(selected):
 _,p,q,h,s,lnN,mode=row;p,q=int(p),int(q)
 A=q*(121*q-99*p)*(26*p-23*q)*(75*p-44*q)
 B=p*(21*q-74*p)*(25*p-36*q)*(p+77*q)
 count=0
 for pr,table in tables.items():
  if q%pr:
   a=table[p*pow(q,-1,pr)%pr]
   if a!=9999:count+=round(128*a*math.log(pr)/pr)
 assert count==round(128*float(s)),(p,q,count,s)
 body.append(f'E=ellminimalmodel(ellinit([0,{A+B},0,{A*B},0]));print("C{i}=",[log(max(abs(E.c4)^3,E.c6^2))/12,sum(j=1,#ps,if(E.disc%ps[j],ellap(E,ps[j])*log(ps[j])/ps[j],0)),elltors(E)[2]]);')
 checks.append({'p':p,'q':q,'height_lower_bound':float(h),'score':float(s)})
out=gp.run_gp('\n'.join(body),timeout=60)
for i,c in enumerate(checks):
 h,s,t=ast.literal_eval(out[f'C{i}']);assert c['height_lower_bound']<=h+1e-10,(c,h)
 assert t==[2,2],t
 c.update(height=h,scalar_gp_score=s)
pathlib.Path('submission/research/x400a-validation.json').write_text(json.dumps(checks,indent=1))
print('VALIDATED',len(checks),'both signs',sum(c['p']<0 for c in checks),'max height gap',max(c['height']-c['height_lower_bound'] for c in checks))
