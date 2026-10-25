import ast,json,pathlib,sys
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[2]))
from ec_verify import gp,verify_text,VerifyConfig
root=pathlib.Path("submission/research")
families=json.loads((root/"constructed-families-strong.json").read_text())
reports=[]
for i,f in enumerate(families,1):
 r,s,v=f["roots"]
 body=f"""read("/home/linux/ec-rank/submission/research/z3-family.gp");
F=z3family({r},{s},{v},{f["lam"]},{f["k"]},101/103);
E=ellminimalmodel(F[1],&vc);P=vector(#F[2],i,ellchangepoint(F[2][i],vc));take=List();
for(i=1,#P,trial=concat(Vec(take),[P[i]]);ev=qfjacobi(ellheightmatrix(E,trial))[1];if(vecmin(ev)>10^-20,listput(take,P[i])));
print("A=",[E.a1,E.a2,E.a3,E.a4,E.a6]);print("N=",#take);
for(i=1,#take,print("X",i,"=",take[i][1]);print("Y",i,"=",take[i][2]));
"""
 out=gp.run_gp(body,timeout=30,precision=57)
 a=ast.literal_eval(out["A"]);n=int(out["N"])
 text="# Section independence check at parameter 101/103; unscored research witness\nE "+" ".join(map(str,a))+"\n"
 text+="".join(f'P {out[f"X{j}"]} {out[f"Y{j}"]}\n' for j in range(1,n+1))
 (root/f"family-{i}-sections.ec").write_text(text)
 result=verify_text(text,VerifyConfig(conductor_enabled=False))
 reports.append(dict(family=f,result=result.to_json()))
 print("CERTIFIED_FAMILY_SECTIONS",i,result.certified_rank,result.torsion_key,"H",result.naive_height,flush=True)
(root/"family-section-certificates.json").write_text(json.dumps(reports,indent=1))
