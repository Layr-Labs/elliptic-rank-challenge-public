#!/usr/bin/env python3
"""Published C2 rank9 and C2xC2 rank6 families, independently chosen parameters.
C2: Elkies--Klagsbrun arXiv:2003.00077v1 §9, Appendix A.
C2xC2: Dujella--Peral arXiv:1712.02082v2 §§2.4--2.5.
No board edits; a full novelty-enabled verifier checks potential score gains.
"""
import argparse,json,math,subprocess,time
from pathlib import Path
from collections import Counter
from concurrent.futures import ThreadPoolExecutor,as_completed
ROOT=Path(__file__).resolve().parents[3]
GP=ROOT/'tools/pari/bin/gp';PY=ROOT/'.venv-bench/bin/python'
C2=r'''
c2(m,t)={
 my(u=2*(m^2-m-1)/(m^2+1),aa,bs,bb,xx);
 aa=(u^8-18*u^6+163*u^4-1152*u^2+4096)*t^4+(3*u^7-35*u^5-120*u^3+1536*u)*t^3+(u^8-13*u^6+32*u^4-152*u^2+1536)*t^2+(u^7+3*u^5-156*u^3+672*u)*t+3*u^6-33*u^4+112*u^2-80;
 bs=[(u^2+u-8)*t-u+2,(u^2-u-8)*t-u-2,(u^2-u-8)*t+u^2+u-10,(u^2+u-8)*t-u^2+u+10,(u^2-7*u+8)*t-u^2+u+2,(u^2+7*u+8)*t+u^2+u-2,(u^2+5*u+8)*t+u^2+3*u+2,(u^2-5*u+8)*t-u^2+3*u-2];
 bb=prod(k=1,8,bs[k]);
 xx=[-bs[1]*bs[2]*bs[3]*bs[6],-bs[1]*bs[2]*bs[4]*bs[5],4*bs[1]*bs[2]*bs[5]*bs[6],bs[1]*bs[3]*bs[4]*bs[6],-bs[1]*bs[3]*bs[4]*bs[7],bs[1]*bs[3]*bs[4]*bs[8],bs[1]*bs[3]*bs[5]*bs[6],-bs[1]*bs[5]*bs[6]*bs[7],-(m-1)^2*bs[1]*bs[2]*bs[3]*bs[8]];
 [ellinit([0,2*aa,0,bb,0]),xx]
};
'''
C22=r'''
c22(v)={
 my(w=2*(v^2-1)/(7*(1+v^2)),aa,bb,xx);
 aa=-2*(-169-12020*w^2+678000*w^4-12680000*w^6+80000000*w^8);
 bb=(-1+10*w)^2*(1+10*w)^2*(-1+20*w^2)*(1+80*w^2)*(-31+400*w^2)*(-41+500*w^2)*(9-200*w^2+2000*w^4);
 xx=[(-1+10*w)^2*(1+10*w)^2*(1+80*w^2)*(-41+500*w^2)/9,(-1+10*w)*(1+10*w)*(1+80*w^2)*(-31+400*w^2)*(-41+500*w^2)/9,(1+80*w^2)*(-11+300*w^2)^2*(-41+500*w^2)/49,(1+80*w^2)*(-31+400*w^2)*(9-200*w^2+2000*w^4),9*(-1+10*w)*(1+10*w)*(-1+20*w^2)*(-41+500*w^2),-9*(-1+10*w)^2*(1+10*w)^2*(-1+20*w^2)*(-41+500*w^2)];
 [ellinit([0,aa,0,bb,0]),xx]
};
'''
def worker(i,params,args,out,targets):
    code=C2+C22+'params=List();\n'+''.join('listput(params,'+json.dumps(p)+');\n' for p in params)
    code+=f'ec_hmax={args.height_max};ec_effort={args.effort};ec_alarm={args.alarm};ec_c2={int(args.family=="c2")};\n'
    code+=r'''{
for(i=1,#params,
 my(p=params[i][1],q=params[i][2],mm=params[i][3],t=p/q,F,E,M,ch,H,P,R,S=0.,tm=getabstime());
 F=if(ec_c2,c2(mm,t),c22(t));E=F[1];
 if(E==[]||E.disc==0,print("ROW ",[p,q,mm,"singular"]);next);
 M=ellminimalmodel(E,&ch);H=log(max(abs(M.c4)^3,M.c6^2))/12.;
 if(H>=ec_hmax,print("ROW ",[p,q,mm,"height",H]);next);
 P=vector(#F[2],k,my(xx=F[2][k],ys=ellordinate(E,xx));if(#ys==0,error("family section invalid"));ellchangepoint([xx,ys[1]],ch));
 setrand(1);R=iferr(alarm(ec_alarm,ellrank(M,ec_effort,P)),e,0);
 if(type(R)!="t_VEC",print("ROW ",[p,q,mm,"timeout",H]);next);
 print("ROW ",[p,q,mm,"rank",H,elltors(M)[2],R[1],R[2],apply(x->Str(x),M[1..5]),apply(P->apply(x->Str(x),P),R[4]),(getabstime()-tm)/1000.]);
);
}
print("COMPLETE");quit;
'''
    cp=subprocess.Popen([str(GP),'-q','-f','--default','parisize=134217728','--default','parisizemax=1073741824','--default','recover=0'],stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,cwd=ROOT)
    assert cp.stdin and cp.stdout and cp.stderr
    cp.stdin.write(code);cp.stdin.close();rows=[];complete=False;hits=[]
    with (out/f'part-{i:05}.jsonl').open('w') as f:
        for line in cp.stdout:
            if line.strip()=='COMPLETE':complete=True
            if not line.startswith('ROW '):continue
            d=json.loads(line[4:]);r=dict(p=d[0],q=d[1],m=d[2],status=d[3],family=args.family)
            if d[3] in ['height','timeout']:r['height']=d[4]
            if d[3]=='rank':
                r.update(height=d[4],torsion=d[5],lower=d[6],upper=d[7],ainvs=d[8],points=d[9],seconds=d[10]);n=len(d[9]);key=('x'.join(map(str,d[5])) or '1')+'/'+str(n);r['cell']=key;target=targets.get(key)
                if n>0 and (target is None or round(d[4],6)<target):
                    file=out/'candidates'/f'{args.family}-m{d[2]}-p{d[0]}-q{d[1]}-r{n}.ec';file.parent.mkdir(exist_ok=True)
                    text='# Independently enumerated '+args.family+' parameter '+str(d[:3])+'\nE '+' '.join(d[8])+'\n'+''.join('P '+' '.join(P)+'\n' for P in d[9])+f'#CLAIM rank={n}\n';file.write_text(text)
                    v=subprocess.run([str(PY),'-m','ec_verify',str(file)],capture_output=True,text=True,timeout=100,cwd=ROOT);file.with_suffix('.verify.txt').write_text(v.stdout+v.stderr);r['verification']={'file':str(file),'exit':v.returncode}
                    if v.returncode==0:
                        z=json.loads(v.stdout);r['verification'].update(cell=z['cell'],height=z['naive_height'],contribution=z['contribution'])
                        if z['contribution']>0:hits.append(r['verification']);print(json.dumps({'stage':'HIT','result':r['verification']}),flush=True)
                    else:r['verification']['rejection']=v.stdout.strip()[:350]
            f.write(json.dumps(r)+'\n');f.flush();rows.append(r)
    stderr=cp.stderr.read();rc=cp.wait();errors=[x for x in stderr.splitlines() if x.lstrip().startswith('***') and 'Warning:' not in x];(out/f'part-{i:05}.stderr.txt').write_text(stderr)
    return dict(index=i,processed=len(rows),expected=len(params),complete=complete and rc==0,errors=errors,hits=hits,ranks=dict(Counter(len(x.get('points',[])) for x in rows if x['status']=='rank')))
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--family',choices=['c2','c22'],required=True);ap.add_argument('--bound',type=int,default=50);ap.add_argument('--m',default='2,-3');ap.add_argument('--height-max',type=float,default=18.998346);ap.add_argument('--workers',type=int,default=20);ap.add_argument('--effort',type=int,default=1);ap.add_argument('--alarm',type=int,default=10);ap.add_argument('--run',required=True);args=ap.parse_args()
    out=Path(args.run);out.mkdir(parents=True,exist_ok=True)
    if list(out.glob('part-*.jsonl')):raise SystemExit('Existing run')
    ms=list(map(int,args.m.split(','))) if args.family=='c2' else [0]
    params=[(p,q,m) for q in range(1,args.bound+1) for p in range(-args.bound,args.bound+1) if p and math.gcd(p,q)==1 for m in ms]
    if args.family=='c22':params=[(p,q,m) for p,q,m in params if p>q]
    targets={k:c['record']['naive_height'] for k,c in json.loads((ROOT/'ec_verify/records.json').read_text())['cells'].items()}
    for k,c in json.loads((ROOT/'score.json').read_text())['metrics']['cells'].items():
        if c.get('board_best') and c.get('contribution',0)>0:targets[k]=c['board_best']['naive_height']
    (out/'parameters.json').write_text(json.dumps(params));(out/'config.json').write_text(json.dumps(vars(args),indent=2));chunks=[params[i:i+16] for i in range(0,len(params),16)];results=[];t0=time.monotonic();print(json.dumps({'stage':'start','parameters':len(params),'batches':len(chunks)}),flush=True)
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        fs=[pool.submit(worker,i,p,args,out,targets) for i,p in enumerate(chunks)]
        for f in as_completed(fs):
            r=f.result();results.append(r)
            if len(results)%10==0 or r['errors'] or r['hits']:print(json.dumps({'stage':'progress','done':len(results),'total':len(chunks),'latest':r}),flush=True)
    rows=[json.loads(s) for p in sorted(out.glob('part-*.jsonl')) for s in p.read_text().splitlines()];(out/'ranked.jsonl').write_text(''.join(json.dumps(r)+'\n' for r in rows));summary=dict(expected=len(params),processed=len(rows),complete=len(rows)==len(params) and all(r['complete'] and not r['errors'] for r in results),states=dict(Counter(r['status'] for r in rows)),ranks=dict(Counter(len(r.get('points',[])) for r in rows if r['status']=='rank')),hits=[h for r in results for h in r['hits']],errors=[r for r in results if r['errors'] or not r['complete']],seconds=time.monotonic()-t0);(out/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps({'stage':'complete',**summary}),flush=True)
if __name__=='__main__':main()
