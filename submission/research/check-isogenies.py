import ast,json,pathlib,sys
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[2]))
from ec_verify import gp,cells,records,verify_text,VerifyConfig
root=pathlib.Path("submission/research")
parent=json.loads((root/"initial-board.json").read_text())
held={r["cell"]:r["contribution"] for r in parent["metrics"]["curves"]}
rows=json.loads((root/"isogeny-upgrade-candidates.json").read_text())
best=[]
for index,row in enumerate(rows):
 c=row["candidate"];pts=row["result"]["points"];n=len(pts)
 pv="["+",".join("["+x+","+y+"]" for x,y in pts)+"]"
 body=f"""E=ellinit({c["a"]});P={pv};L=ellisomat(E,2)[1];print("N=",#L);
for(j=1,#L,Em=ellminimalmodel(ellinit(L[j][1]),&vc);print("A",j,"=",[Em.a1,Em.a2,Em.a3,Em.a4,Em.a6]);print("T",j,"=",elltors(Em)[2]);print("H",j,"=",log(max(abs(Em.c4)^3,Em.c6^2))/12);Q=vector(#P,k,ellchangepoint(ellisogenyapply(L[j][2],P[k]),vc));for(k=1,#Q,print("X",j,"_",k,"=",Q[k][1]);print("Y",j,"_",k,"=",Q[k][2])));
"""
 out=gp.run_gp(body,timeout=30)
 for iso in range(1,int(out["N"])+1):
  h=round(float(out[f"H{iso}"]),6);tor=ast.literal_eval(out[f"T{iso}"]);key=cells.torsion_key(tor)
  gain,r=max((cells.curve_contribution(r,h,records.baseline(records.load(),f"{key}/{r}"))-held.get(f"{key}/{r}",0),r) for r in range(1,n+1))
  print("ISOGENY",c["u"],c["p"],c["q"],iso,key,n,h,"gain",gain,flush=True)
  if gain>0:
   text=f'# Own specialization {c["u"]}, t={c["p"]}/{c["q"]}; rational isogeny {iso}\n'
   a=ast.literal_eval(out[f"A{iso}"]);text+="E "+" ".join(map(str,a))+"\n"
   text+="".join(f'P {out[f"X{iso}_{k}"]} {out[f"Y{iso}_{k}"]}\n' for k in range(1,r+1))
   best.append((gain,c,iso,text))
if best:
 gain,c,iso,text=max(best,key=lambda x:x[0])
 result=verify_text(text,VerifyConfig(conductor_enabled=False))
 name=f'first-z4-isogeny-{c["p"]}-{c["q"]}-{iso}.ec'
 pathlib.Path("submission/board",name).write_text(text,encoding="ascii")
 (root/"discovery.json").write_text(json.dumps(dict(candidate=c,isogeny=iso,result=result.to_json(),delta=gain,file=name),indent=1))
 print("DISCOVERY",name,result.cell,result.naive_height,gain,flush=True)
else:print("NO_ISOGENY_IMPROVEMENT",flush=True)
