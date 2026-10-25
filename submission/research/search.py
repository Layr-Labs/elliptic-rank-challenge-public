"""Independent rational-parameter search in Tate normal form.
Run from the repository root using .venv-bench/bin/python.
Only the official verifier can accept a candidate; heuristics only order work.
"""
import ast
import argparse
import concurrent.futures
import collections
import json
import pathlib
import sys
import time
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))
from ec_verify import VerifyConfig, VerifyError, verify_text
from ec_verify import cells, fingerprint, gp, known_curves, records

AP = argparse.ArgumentParser()
AP.add_argument("--input",default="sieve.txt")
AP.add_argument("--limit",type=int,default=60)
ARGS = AP.parse_args()
ROOT = pathlib.Path("submission/research")
TABLE = records.load()
INITIAL = json.loads((ROOT/"initial-board.json").read_text())
PARENT = {x["cell"]: x["contribution"] for x in INITIAL["metrics"]["curves"]}
import hashlib
parent_digest=hashlib.sha256((ROOT/"initial-board.json").read_bytes()).hexdigest()
cache_path=pathlib.Path("/tmp/ec-rank-completed-curve-checks.json")
try:
 saved_cache=json.loads(cache_path.read_text())
except (OSError,ValueError):
 saved_cache={}
completed=saved_cache.get("checks",{}) if saved_cache.get("parent_digest")==parent_digest else {}
history=ROOT/"search-log.jsonl"
if history.exists():
 for line in history.read_text().splitlines():
  try:
   old=json.loads(line);oldkey=",".join(map(str,old["candidate"]["a"]))
   if old["status"] in ("known","board_duplicate"):
    completed[oldkey]=old["status"]
   elif old["status"]=="rank" and old["result"]["bounds"][1]<old["result"]["target"]:
    completed[oldkey]="upper_bound_below_target"
  except (ValueError,KeyError,TypeError):
   continue

CFG = VerifyConfig(conductor_enabled=False)
GPBIN = gp.find_gp()
start = time.monotonic()
buckets = collections.defaultdict(list)
for line in (ROOT/ARGS.input).read_text().splitlines():
    if not line.startswith("["): continue
    values=ast.literal_eval(line)
    f,p,q,h,s,a=values[:6]
    u=values[6] if len(values)>6 else None
    iso=values[7] if len(values)>7 else 0
    target = next(r for r in range(1,33) if records.baseline(TABLE,f"{f}/{r}") is None
                  or round(h,6) < records.baseline(TABLE,f"{f}/{r}"))
    if ",".join(map(str,a)) in completed:continue
    buckets[f,target].append(dict(f=f,p=p,q=q,h=h,s=s,a=a,target=target,**({"u":u,"iso":iso} if u is not None else {})))
jobs = []
for key, group in buckets.items():
    if key[1] < 2:
        continue
    jobs.extend(sorted(group,key=lambda x:x["s"])[:ARGS.limit])
# Favor rank 3-4 in less represented torsion groups.
jobs.sort(key=lambda x: ((0 if x["f"]==9 else 1 if x["f"]==7 else 2),x["target"]-0.12*(-x["s"]),x["s"]))
(ROOT/"shortlist.json").write_text(json.dumps(jobs,indent=1))
print("SHORTLIST", len(jobs), flush=True)

def investigate(c):
    a = c["a"]
    body = f"""setrand(20261025);
E=ellinit({a});
ps=primes([2,1000]);
print("FP=",vector(#ps,i,if(E.disc%ps[i],abs(ellap(E,ps[i])),-1)));
print("TOR=",elltors(E)[2]);
"""
    try:
        preliminary=gp.run_gp(body,timeout=8,gp_bin=GPBIN)
        fp=ast.literal_eval(preliminary["FP"])
        if known_curves.match(fp) is not None:
            return c, "known", None
        if any(fingerprint.matches(fp, x["fingerprint"]) for x in INITIAL["metrics"]["curves"]):
            return c, "board_duplicate", None
        tors=ast.literal_eval(preliminary["TOR"])
        tkey=cells.torsion_key(tors)
        target=next(r for r in range(1,33) if cells.curve_contribution(r,round(c["h"],6),records.baseline(TABLE,f"{tkey}/{r}"))>PARENT.get(f"{tkey}/{r}",0))
        seed=""
        if "u" in c:
            if c["u"]=="x400a":
                seed=f'read("/home/linux/ec-rank/submission/research/x400a-family.gp");F=x400afamily({c["p"]},{c["q"]});origE=F[1];P0=F[2];'
            elif c["u"].startswith("z4:"):
                family=c["u"][3:]
                seed=f'read("/home/linux/ec-rank/submission/research/z4-family.gp");F=z4family({family},{c["p"]}/{c["q"]});origE=F[1];P0=F[2];'
            elif c["u"].startswith("z3:"):
                params=c["u"][3:]
                seed=f'read("/home/linux/ec-rank/submission/research/z3-family.gp");F=z3family({params},{c["p"]}/{c["q"]});origE=F[1];P0=F[2];'
            else:
                seed=f'read("/home/linux/ec-rank/submission/research/z2-family.gp");F=z2family({c["u"]},{c["p"]}/{c["q"]});origE=F[1];P0=F[2];'
            if c.get("iso"):
                seed+='P0=vector(#P0,i,if(P0[i][1]==0,[0],[(P0[i][2]/P0[i][1])^2,P0[i][2]*(1-origE.a4/P0[i][1]^2)]));origE=ellinit([0,-2*origE.a2,0,origE.a2^2-4*origE.a4,0]);'
            seed+='Em=ellminimalmodel(origE,&v);P=vector(#P0,i,ellchangepoint(P0[i],v));if(!vecprod(vector(#P,i,ellisoncurve(E,P[i]))),error("seed transform failed"));'
        body=f"""setrand(20261025);
E=ellinit({a});
{seed}
R=ellrank(E,0{",P" if seed else ""});
if(R[2]>={target} && #R[4]<{target},R=ellrank(E,4,R[4]));
print("BOUNDS=",[R[1],R[2]]);
print("N=",#R[4]);
for(i=1,#R[4],print("X",i,"=",R[4][i][1]);print("Y",i,"=",R[4][i][2]));
"""
        out=gp.run_gp(body,timeout=45,gp_bin=GPBIN)
        n=int(out["N"])
        pts=[(out[f"X{i}"],out[f"Y{i}"]) for i in range(1,n+1)]
        return c,"rank",dict(bounds=ast.literal_eval(out["BOUNDS"]),n=n,points=pts,tkey=tkey,target=target)
    except (gp.GPError,StopIteration) as e:
        return c,"error",str(e)

counts=collections.Counter()
log=open(ROOT/"search-log.jsonl","a",buffering=1)
found=False
pool=concurrent.futures.ThreadPoolExecutor(max_workers=4)
pending={}
iterator=iter(jobs)
for _ in range(4):
    c=next(iterator,None)
    if c:
        pending[pool.submit(investigate,c)]=c
while pending and not found:
    finished,_=concurrent.futures.wait(pending,return_when=concurrent.futures.FIRST_COMPLETED)
    for task in finished:
        pending.pop(task)
        c,status,r=task.result()
        counts[status]+=1
        cachekey=",".join(map(str,c["a"]))
        if status in ("known","board_duplicate") or (status=="rank" and r["bounds"][1]<r["target"]):
            completed[cachekey]=status
        if sum(counts.values())%100==0:
            cache_path.write_text(json.dumps(dict(parent_digest=parent_digest,checks=completed)))

        log.write(json.dumps(dict(candidate=c,status=status,result=r))+"\n")
        if status=="rank":
            print("RANK",c["f"],f'{c["p"]}/{c["q"]}',"H",round(c["h"],6),r["bounds"],"points",r["n"],"target",r["target"],flush=True)
            options=[]
            for n in range(1,r["n"]+1):
                cell=f'{r["tkey"]}/{n}'
                value=cells.curve_contribution(n,round(c["h"],6),records.baseline(TABLE,cell))
                options.append((value-PARENT.get(cell,0),n))
            if options and max(options)[0]>0:
                delta,n=max(options)
                text=f'# Found by independent Tate-family search; family={c["f"]}; t={c["p"]}/{c["q"]}\n'
                text+=f'# family u={c.get("u","not-applicable")}\n'
                text+="E "+" ".join(map(str,c["a"]))+"\n"
                text+="".join(f"P {x} {y}\n" for x,y in r["points"][:n])
                try:
                    result=verify_text(text,CFG)
                    filename=f'first-t{c["f"]}-p{c["p"]}-q{c["q"]}-r{n}.ec'
                    pathlib.Path("submission/board",filename).write_text(text,encoding="ascii")
                    (ROOT/"discovery.json").write_text(json.dumps(dict(candidate=c,result=result.to_json(),delta=delta,file=filename),indent=1))
                    print("DISCOVERY",filename,"CELL",result.cell,"H",result.naive_height,"GAIN",delta,flush=True)
                    found=True
                except VerifyError as e:
                    print("VERIFICATION_REJECTED",str(e),flush=True)
        if sum(counts.values())%25==0:
            print("PROGRESS",dict(counts),"seconds",round(time.monotonic()-start,1),flush=True)
        if not found:
            c=next(iterator,None)
            if c:
                pending[pool.submit(investigate,c)]=c
for task in pending:
    task.cancel()
pool.shutdown(wait=True,cancel_futures=True)
cache_path.write_text(json.dumps(dict(parent_digest=parent_digest,checks=completed)))
print("FINISHED",dict(counts),"found",found,"seconds",round(time.monotonic()-start,1),flush=True)
