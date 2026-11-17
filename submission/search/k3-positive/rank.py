#!/usr/bin/env python3
"""Exact seeded ranking of our own K3 specializations; full verifier on hits.
Family and four sections: Elkies, as reproduced in arXiv:2207.08206v2 §3.
Sieve parameter pairs must be generated independently. Never edits board.
"""
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
from pathlib import Path
import subprocess
import time

ROOT=Path(__file__).resolve().parents[3]
GP=ROOT/'tools/pari/bin/gp'
PY=ROOT/'.venv-bench/bin/python'

def batch(index, pairs, args, run, baselines):
    code='params=List();\n'+''.join('listput(params,'+json.dumps(pair)+');\n' for pair in pairs)
    code+=f'ec_alarm={args.alarm};ec_effort={args.effort};ec_hmax={args.height_max};ec_second={int(args.second)};\n'
    code+=r'''{
for(i=1,#params,
 my(p=params[i][1],q=params[i][2],t=p/q,aa=(8*t-1)*(32*t+7),bb=8*(t+1)*(15*t-8)*(31*t-7),E,M,ch,xx,pts,R,H,tm=getabstime());
 if(ec_second,aa=4*(9+80*t);bb=(t-2)*(2*t-81)*(2*t-1)*(18*t-1)/2);
 E=ellinit([0,aa^2-8*bb,0,16*bb^2,0]);
 if(E==[]||E.disc==0,print("ROW ",[p,q,"singular"]);next);
 M=ellminimalmodel(E,&ch);H=log(max(abs(M.c4)^3,M.c6^2))/12.;
 if(H>=ec_hmax,print("ROW ",[p,q,"height",H]);next);
 if(ec_second,pts=[],
   xx=[-361*(t+1)*(31*t-7),-4*(t+1)*(15*t-8)*(16*t-7)^2,-16*(t+1)*(8*t+7)^2*(15*t-8),4*(15*t-8)*(16*t+1)^2*(31*t-7)];
   pts=vector(4,k,my(ys=ellordinate(E,xx[k]));if(#ys==0,error("section off curve"));ellchangepoint([xx[k],ys[1]],ch));
 );
 setrand(1);R=iferr(alarm(ec_alarm,ellrank(M,ec_effort,pts)),e,0);
 if(type(R)!="t_VEC",print("ROW ",[p,q,"timeout",H]);next);
 print("ROW ",[p,q,"rank",H,elltors(M)[2],R[1],R[2],apply(x->Str(x),M[1..5]),apply(P->apply(x->Str(x),P),R[4]),(getabstime()-tm)/1000.]);
);
}
print("COMPLETE");quit;
'''
    proc=subprocess.Popen([str(GP),'-q','-f','--default','parisize=134217728','--default','parisizemax=1073741824','--default','recover=0'],stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,cwd=ROOT)
    assert proc.stdin is not None and proc.stdout is not None and proc.stderr is not None
    proc.stdin.write(code);proc.stdin.close();complete=False;rows=[]
    with (run/f'part-{index:05}.jsonl').open('w') as out:
        for line in proc.stdout:
            if line.strip()=='COMPLETE':complete=True
            if not line.startswith('ROW '):continue
            v=json.loads(line[4:]);r={'p':v[0],'q':v[1],'status':v[2]}
            if v[2] in ['height','timeout']:r['height']=v[3]
            if v[2]=='rank':
                r.update(height=v[3],torsion=v[4],lower=v[5],upper=v[6],ainvs=v[7],points=v[8],seconds=v[9])
                tor='x'.join(map(str,r['torsion'])) or '1';n=len(r['points']);cell=f'{tor}/{n}';r['cell']=cell
                best=baselines.get(cell)
                if n>0 and (best is None or round(r['height'],6)<best):
                    p=run/'candidates'/f'k3-p{r["p"]}-q{r["q"]}-r{n}.ec';p.parent.mkdir(exist_ok=True)
                    text='# Original K3 specialization '+str((r['p'],r['q']))+'\nE '+' '.join(r['ainvs'])+'\n'+''.join('P '+' '.join(P)+'\n' for P in r['points'])+f'#CLAIM rank={n}\n';p.write_text(text)
                    cp=subprocess.run([str(PY),'-m','ec_verify',str(p)],cwd=ROOT,capture_output=True,text=True,timeout=90)
                    p.with_suffix('.verify.txt').write_text(cp.stdout+cp.stderr)
                    r['verification']={'file':str(p),'exit':cp.returncode}
                    if cp.returncode==0:
                        d=json.loads(cp.stdout);r['verification'].update(cell=d['cell'],height=d['naive_height'],contribution=d['contribution'])
                        print(json.dumps({'stage':'VERIFIED','result':r['verification']}),flush=True)
                    else:r['verification']['rejection']=(cp.stdout+cp.stderr).strip()[:450]
            out.write(json.dumps(r)+'\n');out.flush();rows.append(r)
    err=proc.stderr.read();rc=proc.wait();errors=[s for s in err.splitlines() if s.lstrip().startswith('***') and 'Warning:' not in s]
    (run/f'part-{index:05}.stderr.txt').write_text(err)
    return {'index':index,'count':len(rows),'expected':len(pairs),'complete':complete and rc==0,'errors':errors,'rank_counts':dict(Counter(len(r.get('points',[])) for r in rows if r['status']=='rank')),'hits':[r['verification'] for r in rows if r.get('verification',{}).get('contribution',0)>0]}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input',required=True);ap.add_argument('--run',required=True);ap.add_argument('--workers',type=int,default=20);ap.add_argument('--batch-size',type=int,default=32);ap.add_argument('--alarm',type=int,default=10);ap.add_argument('--effort',type=int,default=1);ap.add_argument('--height-max',type=float,default=18.746574);ap.add_argument('--retry',action='store_true');ap.add_argument('--second',action='store_true');args=ap.parse_args()
    run=Path(args.run);run.mkdir(parents=True,exist_ok=True)
    if list(run.glob('part-*.jsonl')):raise SystemExit('Refusing to overwrite run')
    if args.retry:
        old=[json.loads(s) for s in Path(args.input).read_text().splitlines()]
        pairs=sorted({(r['p'],r['q']) for r in old if r['status']=='timeout' or (r['status']=='rank' and r.get('upper',0)>len(r.get('points',[])))})
    else:pairs=sorted({tuple(map(int,s.split()[:2])) for s in Path(args.input).read_text().splitlines() if s.strip()})
    records=json.loads((ROOT/'ec_verify/records.json').read_text())['cells'];baselines={k:v['record']['naive_height'] for k,v in records.items()}
    score=json.loads((ROOT/'score.json').read_text())
    for k,c in score['metrics']['cells'].items():
        if c.get('board_best') and c.get('contribution',0)>0:baselines[k]=c['board_best']['naive_height']
    (run/'config.json').write_text(json.dumps(vars(args),indent=2));(run/'parameters.json').write_text(json.dumps(pairs));t0=time.monotonic();chunks=[pairs[i:i+args.batch_size] for i in range(0,len(pairs),args.batch_size)];results=[]
    print(json.dumps({'stage':'start','pairs':len(pairs),'chunks':len(chunks)}),flush=True)
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        fs=[pool.submit(batch,i,p,args,run,baselines) for i,p in enumerate(chunks)]
        for f in as_completed(fs):
            result=f.result();results.append(result)
            if len(results)%10==0 or result['hits'] or result['errors']:print(json.dumps({'stage':'progress','done':len(results),'total':len(chunks),'latest':result}),flush=True)
    rows=[json.loads(s) for p in sorted(run.glob('part-*.jsonl')) for s in p.read_text().splitlines()];(run/'ranked.jsonl').write_text(''.join(json.dumps(r)+'\n' for r in rows))
    summary={'expected':len(pairs),'processed':len(rows),'complete':len(rows)==len(pairs) and all(r['complete'] and not r['errors'] and r['count']==r['expected'] for r in results),'states':dict(Counter(r['status'] for r in rows)),'rank_counts':dict(Counter(len(r.get('points',[])) for r in rows if r['status']=='rank')),'hits':[h for r in results for h in r['hits']],'errors':[r for r in results if r['errors'] or not r['complete']],'seconds':time.monotonic()-t0};(run/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps({'stage':'complete',**summary}),flush=True)

if __name__=='__main__':main()
