import struct,sys,subprocess,collections
def vl(x,l):
    k=0
    while x%l==0 and k<2: x//=l;k+=1
    return k
def M_of(b): return (2 if vl(b,2)==1 else 1)*(3 if vl(b,3)==1 else 1)
def rows(f):
    d=open(f,'rb').read(); return [struct.unpack_from('<qqii',d,24*k) for k in range(len(d)//24)]
def spans(f):
    d=open(f,'rb').read(); return [struct.unpack_from('<iii',d,12*k) for k in range(len(d)//12)]
SR=rows('o_o16.rows64'); SS=spans('o_o16.spans')
bv=lambda r:(r[1]<<64)|(r[0]&((1<<64)-1))
for si,so in [(0,0),(1234,0),(1234,777777),(len(SS)-1,5),(2000,123456789)]:
    if si>=len(SS): continue
    so=min(so,SS[si][2]-1)
    open('ck.txt','w').write(f"{si} {so} 1000 7 12.5\n")
    out=subprocess.run(['../port/vslice','v_o16','vs_test','o_o16','ck.txt'],capture_output=True,text=True)
    r0=SS[si][0]; b_r=bv(SR[r0]); a_start=SS[si][1]+so
    exp=collections.defaultdict(list)
    for j in range(si,len(SS)):
        r,i0,ln=SS[j]; lo,hi=i0,i0+ln-1
        if j==si: lo=a_start
        b=bv(SR[r]); M=M_of(b); c0=-((-lo)//M); c1=hi//M
        if c1>=c0: exp[b].append((c0,c1))
    VR=rows('vs_test.vrows64'); VS=spans('vs_test.spans'); got=collections.defaultdict(list)
    for r,i0,ln in VS: got[bv(VR[r])].append((i0,i0+ln-1))
    def norm(d):
        o={}
        for b,L in d.items():
            L=sorted(L); m=[]
            for a,c in L:
                if m and a<=m[-1][1]+1: m[-1]=(m[-1][0],max(m[-1][1],c))
                else: m.append((a,c))
            o[b]=m
        return o
    ok=norm(exp)==norm(got) and all(VR[k][2]==M_of(bv(VR[k])) for k in range(len(VR)))
    print(f"ck si={si} so={so}: vslice says '{out.stdout.strip()}' ; rows {len(VR)} ; match {ok}")
