# check: lotenum/lotenumwd --vmul spans == unrestricted spans filtered by (v2(b)=1 -> 2|a, v3(b)=1 -> 3|a)
import struct, sys, collections
def vl(x,l):
    k=0
    while x%l==0 and k<2: x//=l; k+=1
    return k
def rule(b): return (2 if vl(b,2)==1 else 1)*(3 if vl(b,3)==1 else 1)
def load(prefix, rowsext, wide):
    rb=open(prefix+rowsext,'rb').read(); sb=open(prefix+'.spans','rb').read()
    rs=24 if wide else 16; rows=[]
    for k in range(len(rb)//rs):
        if wide: lo,hi,tri,pad=struct.unpack_from('<qqii',rb,rs*k); b=(hi<<64)|(lo & ((1<<64)-1)); rows.append((b,tri,pad))
        else: q,M,tri=struct.unpack_from('<qii',rb,rs*k); rows.append((q,M,tri))
    iv=collections.defaultdict(list)
    for k in range(len(sb)//12):
        r,i0,ln=struct.unpack_from('<iii',sb,12*k); iv[r].append((i0,i0+ln-1))
    out={}
    for r,(b,M,x) in enumerate(rows):
        L=sorted(iv[r]); m=[]
        for a,c in L:
            if m and a<=m[-1][1]+1: m[-1]=(m[-1][0],max(m[-1][1],c))
            else: m.append((a,c))
        out[b]=(M,x,m)
    return out
name,wide=sys.argv[1],sys.argv[2]=='w'
O=load('o_'+name,'.rows64' if wide else '.rows',wide); V=load('v_'+name,'.vrows64' if wide else '.vrows',wide)
bad=0; nrow=0; kept=0; full=0; Mc=collections.Counter()
for b,(M0,x0,iv0) in O.items():
    M=rule(b); Mc[M]+=1
    if not wide: assert M0==1
    exp=[]
    for lo,hi in iv0:
        full+=hi-lo+1
        i0=-((-lo)//M); i1=hi//M
        if i1>=i0: exp.append((i0,i1))
    got=V.get(b)
    if got is None:
        if exp: bad+=1; print("missing row",b)
        continue
    gM=got[0] if not wide else got[0]   # z3b: M field; z3w: tri
    if gM!=M or (wide and got[1]!=0) or got[2]!=exp: bad+=1; print("MISMATCH",b,M,got[:2],got[2][:3],exp[:3]) if bad<5 else None
    nrow+=1; kept+=sum(c-a+1 for a,c in got[2])
extra=[b for b in V if b not in O]
print(f"{name}: rows {len(O)} vrows {len(V)} extra {len(extra)} mismatches {bad}  M counts {dict(Mc)}  kept {kept} full {full} tau {kept/full:.6f}")
# brute force on the 60 rows with the smallest a-range: expand pairs
NB=int(sys.argv[3]) if len(sys.argv)>3 else 60
small=sorted(O.items(), key=lambda kv: sum(c-a+1 for a,c in kv[1][2]))[:NB]; bf=0
for b,(M0,x0,iv0) in small:
    M=rule(b); A=set()
    for lo,hi in iv0:
        for a in range(lo,hi+1):
            if (vl(b,2)!=1 or a%2==0) and (vl(b,3)!=1 or a%3==0): A.add(a)
    G=set(); g=V.get(b)
    if g:
        for i0,i1 in g[2]:
            for i in range(i0,i1+1): G.add(g[0]*i)
    if A!=G: bf+=1
print(f"  brute force: {len(small)} rows, {sum(1 for _ in small)} compared, set mismatches {bf}")
