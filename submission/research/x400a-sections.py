import itertools,json,pathlib
import sympy as S
t,c=S.symbols('t c')
aa=[S.Integer(1),121-99*t,26*t-23,75*t-44]
bb=[t,21-74*t,25*t-36,t+77]
sections=set()
for ai in itertools.combinations(range(4),2):
 for bi in itertools.combinations(range(4),2):
  al=S.prod(aa[i] for i in ai);ar=S.prod(aa[i] for i in range(4) if i not in ai)
  bl=S.prod(bb[i] for i in bi);br=S.prod(bb[i] for i in range(4) if i not in bi)
  f=S.Poly(c*bl+ar,t);g=S.Poly(c*al+br,t)
  df=S.discriminant(f.as_expr(),t);dg=S.discriminant(g.as_expr(),t)
  common=S.gcd(df,dg)
  roots=S.solve(common,c) if S.degree(common,c)>0 else []
  for val in roots:
   if not val.is_Rational or not val:continue
   x=S.expand(val*al*bl);v=S.factor(x*(x+S.prod(aa))*(x+S.prod(bb)))
   y=S.sqrtdenest(S.sqrt(v))
   # Direct squarefree factor check for a rational polynomial y.
   coef,factors=S.factor_list(v,t)
   if any(e%2 for _,e in factors) or not S.sqrt(coef).is_Rational:continue
   y=S.sqrt(coef)*S.prod(k**(e//2) for k,e in factors)
   assert S.expand(y*y-v)==0
   sections.add((str(S.factor(x)),str(S.factor(y))))
result=sorted(sections)
pathlib.Path('submission/research/x400a-sections.json').write_text(json.dumps(result,indent=1))
print('SECTIONS',len(result))
for x,y in result:print(x,';',y)
