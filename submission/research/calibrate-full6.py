"""Known baseline is a control ONLY, never a submission candidate."""
import json,pathlib,sys,ast
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[2]))
from ec_verify import gp
root=pathlib.Path("submission/research")
rec=json.load(open("ec_verify/records.json"))["cells"]["6x2/6"]["record"]
a=[int(x) for x in rec["ainvs"]]
body=f"""E0=ellinit({a});J0=E0.j;
A=37-84*x+102*x^2-36*x^3-3*x^4;B=32*(x-1)^3*(x+1)^3*(3*x-5);c4=16*(A^2-3*B);D=16*B^2*(A^2-4*B);pol=numerator(c4^3-J0*D);fac=factor(pol);roots=List();
for(i=1,matsize(fac)[1],ff=fac[i,1];if(poldegree(ff)==1,listput(roots,-polcoeff(ff,0)/polcoeff(ff,1))));
print("N=",#roots);for(i=1,#roots,print("P",i,"=",numerator(roots[i]));print("Q",i,"=",denominator(roots[i])));
"""
out=gp.run_gp(body,timeout=60)
params=[(int(out[f"P{i}"]),int(out[f"Q{i}"])) for i in range(1,int(out["N"])+1)]
params.sort(key=lambda x:max(abs(x[0]),x[1]))
report=dict(known_control_id=rec["id"],baseline_height=rec["naive_height"],parameters=params,not_a_submission=True)
(root/"full6-baseline-control.json").write_text(json.dumps(report,indent=1))
print("KNOWN_BASELINE_PARAMETER_SIZES",[max(abs(p),q) for p,q in params],flush=True)
