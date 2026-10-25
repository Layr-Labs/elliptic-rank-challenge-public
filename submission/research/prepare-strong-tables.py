import json,pathlib
from fractions import Fraction
root=pathlib.Path("submission/research")
families=json.loads((root/"constructed-families-strong.json").read_text())
body=['{']
rows=[]
for i,f in enumerate(families,1):
 r,s,v=f["roots"];m=f["m"];b=int(m*Fraction(f["lam"]));k=int(m**3*Fraction(f["k"]))
 rows.append(f"{i} {r} {s} {v} {m} {b} {k}")
 body += [f'forprime(p=2,1000,vap=vector(p);for(j=0,p-1,t=Mod(j,p);E=ellinit([{m}*t^2+({b}),0,({k})*(t^2-{r*r})*(t^2-{s*s})*(t^2-{v*v}),0,0]);vap[j+1]=if(#E,p+1-ellcard(E),9999));Einf=ellinit([Mod({m},p),0,Mod({k},p),0,0]);apinf=if(#Einf,p+1-ellcard(Einf),9999);print({i}," ",p," ",concat(vap,[apinf])));']
body+=['}','quit;']
(root/"strong-coefficients.txt").write_text("\n".join(rows)+"\n")
(root/"strong-tables.gp").write_text("\n".join(body)+"\n")
