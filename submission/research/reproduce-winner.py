"""Regenerate the winning points from a fresh family parameter, then verify."""
import ast,json,pathlib,sys
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[2]))
from ec_verify import gp,verify_text,VerifyConfig
root=pathlib.Path(__file__).resolve().parent
family=(root/'x400a-family.gp').as_posix()
out=gp.run_gp(f"""setrand(20261025);
read("{family}");
F=x400afamily(35213,16986);
E=ellminimalmodel(F[1],&v);
P=ellchangepoint(F[2],v);
R=ellrank(E,0,P);
if(R[2]>=13 && #R[4]<13,R=ellrank(E,4,R[4]));
print("A=",[E.a1,E.a2,E.a3,E.a4,E.a6]);
print("BOUNDS=",[R[1],R[2]]);
print("N=",#R[4]);
for(i=1,#R[4],print("X",i,"=",R[4][i][1]);print("Y",i,"=",R[4][i][2]));
""",timeout=120)
assert int(out['N'])>=13,out
a=ast.literal_eval(out['A'])
text='E '+' '.join(map(str,a))+'\n'+''.join(f'P {out[f"X{i}"]} {out[f"Y{i}"]}\n' for i in range(1,14))
board=root.parent/'board'/'first-t2x2-p35213-q16986-r13.ec'
actual='\n'.join(s for s in board.read_text().splitlines() if not s.startswith('#'))+'\n'
assert text==actual,'Fresh seeded point search differs from submitted witnesses'
result=verify_text(text,VerifyConfig(conductor_enabled=False))
(root/'winner-reproduced.json').write_text(json.dumps(dict(bounds=ast.literal_eval(out['BOUNDS']),result=result.to_json()),indent=1))
print('Fresh reproduction:',out['BOUNDS'],result.cell,'H',result.naive_height,'13 identical witnesses; official verification passed')
