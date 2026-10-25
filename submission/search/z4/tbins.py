import math,sys
edges=[-1e30,-30,-10,-3,-1,-0.3,-0.0625,-0.02,0,0.02,0.0625,0.3,1,3,10,1e30]
def bin(t):
    for i in range(len(edges)-1):
        if edges[i]<=t<edges[i+1]: return i
def H4(s,u):
    b2=s*s-4*u; b4=-s*s*u; b6=s*s*u*u
    c4=b2*b2-24*b4; c6=-b2**3+36*b2*b4-216*b6
    return max(3*math.log(abs(c4)+1e-300),2*math.log(abs(c6)+1e-300))/12
lin=lambda a,b,n:[a+(b-a)*i/(n-1) for i in range(n)]
X=float(sys.argv[1]); XL=float(sys.argv[2])
vol=[0]*15
for ls in lin(0,X+1.5,240):
    s=math.exp(ls)
    for sg in (1,-1):
        for lt in lin(-10,25,2400):
            t=sg*math.exp(lt); u=t*s*s; h=H4(s,u)
            if XL<h<X: vol[bin(t)]+=abs(u)*s
tv=sum(vol)
cnt={7:[0]*15,8:[0]*15,6:[0]*15}
for f in sys.argv[3:]:
    for l in open(f):
        x=l.split(); r=int(x[6]); h=float(x[5])
        if r in cnt and XL<h<X: cnt[r][bin(int(x[0])/(int(x[1])*int(x[2])**2))]+=1
for i in range(15):
    print("[%g,%g)"%(edges[i],edges[i+1]), "vol %.3f"%(vol[i]/tv), "r6 %d r7 %d r8 %d"%(cnt[6][i],cnt[7][i],cnt[8][i]), "enr7 %.2f"%(cnt[7][i]/max(1,sum(cnt[7]))/(vol[i]/tv+1e-9)))
