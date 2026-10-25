// X400a Z/2xZ/2 family sieve: y^2 = x(x+A)(x+B),
// A = q(121q-99p)(26p-23q)(75p-44q), B = p(21q-74p)(25p-36q)(p+77q); t=p/q, q in[1,N] (1+k mod K), p in [-N,N]
// usage: x4 N thr1 thr2 Hmax K k [qmin]   output: p q S1 S2 Hest(min over E and 3 2-isogenous) H0est
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
typedef __int128 i128;
static int primes[1000], np=0, np1=0; static float *W[1000];
static long md(long a,long p){a%=p; return a<0?a+p:a;}
static long pw(long a,long e,long p){long r=1;a=md(a,p);while(e){if(e&1)r=r*a%p;a=a*a%p;e>>=1;}return r;}
static long inv(long a,long p){return pw(a,p-2,p);}
static const long FA[4][2]={{0,1},{-99,121},{26,-23},{75,-44}}; // A factors: c_p p + c_q q
static const long FB[4][2]={{1,0},{-74,21},{25,-36},{1,77}};
static const long FD[4][2]={{1,-23},{1,-4},{25,-11},{74,-121}}; // A-B
static void build(int lo,int hi){
  for(int p=lo;p<=hi;p++){int ok=1;for(int d=2;d*d<=p;d++)if(p%d==0)ok=0; if(!ok)continue;
    primes[np]=p; W[np]=calloc(p+1,sizeof(float));
    signed char*chi=malloc(p); memset(chi,-1,p); chi[0]=0; for(long y=1;y<p;y++) chi[y*y%p]=1;
    for(long t=0;t<=p;t++){ long P=(t<p)?t:1, Q=(t<p)?1:0;
      long A=1,B=1; for(int i=0;i<4;i++){ A=A*md(FA[i][0]*P+FA[i][1]*Q,p)%p; B=B*md(FB[i][0]*P+FB[i][1]*Q,p)%p; }
      long s=0; for(long x=0;x<p;x++){ long f=x*((x+A)%p)%p*((x+B)%p)%p; s+=chi[f]; }
      W[np][t]=(float)(s*log((double)p)/p); }
    free(chi); np++; }
}
static int vl(i128 x,long l){ if(x==0) return 99; int v=0; if(x<0)x=-x; while(x%l==0){x/=l;v++;} return v; }
static long double Hraw(long double a2,long double a4){ long double c4=16*(a2*a2-3*a4), c6=-32*a2*(2*a2*a2-9*a4);
  return fmaxl(3*logl(fabsl(c4)),2*logl(fabsl(c6)))/12; }
static double H0g;
static int vm(i128 x, i128 m, long l, int e){ x%=m; if(x<0) x+=m; if(x==0) return e; int v=0; while(x%l==0){x/=l;v++;} return v; }
static double hest(long p,long q){
  i128 A=1,B=1,D=1; for(int i=0;i<4;i++){ A*=(i128)(FA[i][0]*p+FA[i][1]*q); B*=(i128)(FB[i][0]*p+FB[i][1]*q); D*=(i128)(FD[i][0]*p+FD[i][1]*q); }
  long double Al=(long double)A, Bl=(long double)B;
  long double h[4]; h[0]=Hraw(Al+Bl,Al*Bl); h[1]=Hraw(-2*(Al+Bl),(Al-Bl)*(Al-Bl)); h[2]=Hraw(2*(2*Al-Bl),Bl*Bl); h[3]=Hraw(2*(2*Bl-Al),Al*Al);
  static const int SP[]={2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113};
  double r[4]={0,0,0,0};
  for(unsigned k=0;k<sizeof(SP)/sizeof(int);k++){ long l=SP[k]; int vA=vl(A,l), vB=vl(B,l), vD=vl(D,l); int v2=(l==2);
    if(l>3 && vA+vB+vD<3) continue;
    i128 m=1; int e=0; while(m < ((i128)1<<40)){ m*=l; e++; }
    i128 Am=A%m, Bm=B%m; if(Am<0)Am+=m; if(Bm<0)Bm+=m;
    i128 S=(Am+Bm)%m, Dm=(Am-Bm+m)%m, T1=(2*Am-Bm+m)%m, T2=(2*Bm-Am+m)%m;
    // curve 0: c4=16(A^2-AB+B^2), c6=-32(A+B)(2A-B)(A-2B), vDelta=4v2+2(vA+vB+vD)
    int c4v,c6v,dv,kk; double lg=log((double)l);
    #define RED(idx,C4,C6,DV) { c4v=vm(C4,m,l,e); c6v=vm(C6,m,l,e); dv=DV; kk=dv/12; if(c4v/4<kk) kk=c4v/4; if(c6v/6<kk) kk=c6v/6; if(kk>0) r[idx]+=kk*lg; }
    RED(0, 16*((Am*Am%m - Am*Bm%m + Bm*Bm%m + m)%m), 32*(S*T1%m)%m*((Am-2*Bm+2*m)%m), 4*v2+2*(vA+vB+vD));
    // curve 1: a2=-2S, a4=D^2: c4=16(4S^2-3D^2), c6=64 S (8S^2-9D^2), Delta=256 A B D^4
    RED(1, 16*((4*S%m*S%m - 3*Dm%m*Dm%m + 2*m)%m), 64*S%m*((8*S%m*S%m - 9*Dm%m*Dm%m + 2*m)%m), 8*v2+vA+vB+4*vD);
    // curve 2: a2=2T1, a4=B^2: c4=16(4T1^2-3B^2), c6=-64 T1(8T1^2-9B^2), Delta=256 B^4 A D
    RED(2, 16*((4*T1%m*T1%m - 3*Bm%m*Bm%m + 2*m)%m), 64*T1%m*((8*T1%m*T1%m - 9*Bm%m*Bm%m + 2*m)%m), 8*v2+4*vB+vA+vD);
    RED(3, 16*((4*T2%m*T2%m - 3*Am%m*Am%m + 2*m)%m), 64*T2%m*((8*T2%m*T2%m - 9*Am%m*Am%m + 2*m)%m), 8*v2+4*vA+vB+vD);
  }
  double best=1e9; for(int i=0;i<4;i++){ double x=(double)h[i]-r[i]; if(x<best) best=x; }
  H0g=(double)h[0]-r[0];
  return best;
}
int main(int argc,char**argv){
  long N=atol(argv[1]); double thr1=atof(argv[2]), thr2=atof(argv[3]), Hmax=atof(argv[4]); int K=atoi(argv[5]), k=atoi(argv[6]);
  long qmin = argc>7?atol(argv[7]):1;
  build(3,500); np1=np; build(501,3000);
  long *idx=malloc(np*sizeof(long)), *step=malloc(np*sizeof(long)); float base;
  long cnt=0,out=0;
  for(long q=qmin+k;q<=N;q+=K){
    base=0;
    for(int i=0;i<np;i++){ long pr=primes[i]; if(q%pr==0){ step[i]=-1; if(i<np1) base+=W[i][pr]; continue;} step[i]=inv(q%pr,pr); idx[i]=md(-N,pr)*step[i]%pr; }
    for(long p=-N;p<=N;p++){
      float s=base; for(int i=0;i<np1;i++){ if(step[i]>=0){ s+=W[i][idx[i]]; idx[i]+=step[i]; if(idx[i]>=primes[i]) idx[i]-=primes[i]; } }
      if(s<thr1) continue;
      long g=p<0?-p:p, bb=q; while(bb){long r=g%bb;g=bb;bb=r;} if(g!=1) continue;
      double s2=s; for(int i=np1;i<np;i++){ long pr=primes[i]; s2+= (step[i]>=0)? W[i][md(p,pr)*step[i]%pr] : W[i][pr]; }
      cnt++;
      if(s2<thr2) continue;
      double h=hest(p,q); if(h-0.3>Hmax) continue;
      printf("%ld %ld %.3f %.3f %.4f %.4f\n",p,q,s,s2,h,H0g); out++;
    }
  }
  fprintf(stderr,"done cnt=%ld out=%ld\n",cnt,out);
  return 0;
}
