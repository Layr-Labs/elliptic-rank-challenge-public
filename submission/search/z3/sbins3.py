import math,sys
edges=[-1e30,-100,-10,-3,-1,-0.3,-0.1,0,0.02,0.037,0.06,0.1,0.3,1,3,10,100,1e30]
def bn(t):
    for i in range(len(edges)-1):
        if edges[i]<=t<edges[i+1]: return i
def H3(a,b):
    c4=a**4-24*a*b; c6=-a**6+36*a**3*b-216*b*b
    return max(3*math.log(abs(c4)+1e-300),2*math.log(abs(c6)+1e-300))/12
lin=lambda a,b,n:[a+(b-a)*i/(n-1) for i in range(n)]
X=float(sys.argv[1]); XL=float(sys.argv[2])
vol=[0]*(len(edges)-1)
for la in lin(0,X+0.5,240):
    a=math.exp(la)
    for sg in (1,-1):
        for lt in lin(-12,20,2400):
            t=sg*math.exp(lt); b=t*a**3; h=H3(a,b)
            if XL<h<X: vol[bn(t)]+=abs(b)*a
tv=sum(vol)
cnt={}
for f in sys.argv[3:]:
    for l in open(f):
        x=l.split(); r=int(x[5]); h=float(x[4])
        if XL<h<X: cnt.setdefault(r,[0]*len(vol)); cnt[r][bn(int(x[1])/int(x[0])**3)]+=1
for i in range(len(vol)):
    tot7=sum(cnt.get(7,[1])); 
    print("[%g,%g) vol %.3f"%(edges[i],edges[i+1],vol[i]/tv), " ".join("r%d %d"%(r,cnt[r][i]) for r in sorted(cnt) if r>=5), " enr7 %.2f"%(cnt.get(7,[0]*len(vol))[i]/max(1,tot7)/(vol[i]/tv+1e-9)))
