// Mestre ten-root sieve over pair-structured root sets.
// roots: pairs (L_i, L_i+g_i), L_1=0, L_{i+1} > L_i+g_i, L_5+g_5 <= M, g_i<=G
// usage: ms M G Hmax Smin [L2lo L2hi]
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
typedef __int128 i128;
static int PR[512], NPR=0, NPR1=0; static signed char *CHI[512]; static double W[512]; static int MODE=0; static double T2=-1e9;
static int M,G; static double HMAX,SMIN;
static double SCLEAN=0;
static long long mm(i128 z,long long m){ long long r=(long long)(z%m); return r<0?r+m:r; }
static double lscale(i128 r4,i128 r3,i128 r2,i128 r1,i128 r0){
  static const int L[]={2,3,5,7,11,13,17,19,23,29,31,37,41,43,47};
  double ls=0;
  for(int t=0;t<15;t++){ int l=L[t]; int e=0;
    while(1){ int ee=e+1; double lg=6.0*ee*log((double)l); if(lg>42.0) break;
      long long m4=1,m6=1; for(int k=0;k<4*ee;k++) m4*=l; for(int k=0;k<6*ee;k++) m6*=l;
      // need c4 mod m4 and c6 mod m6; compute with modulus m6 (multiple of m4)
      long long m=m6; i128 A=mm(r4,m),B=mm(r3,m),C=mm(r2,m),D=mm(r1,m),E=mm(r0,m);
      i128 I=( (12*A%m)*E%m - (3*B%m)*D%m + C*C%m )%m; I=(I%m+m)%m;
      i128 J=( (((72*A%m)*C%m)*E)%m + (((9*B%m)*C%m)*D)%m - (((27*A%m)*D%m)*D)%m - (((27*E%m)*B%m)*B)%m - (((2*C%m)*C%m)*C)%m )%m; J=(J%m+m)%m;
      i128 c4=(1296*I)%m, c6=(23328*J)%m;
      if(c4%m4!=0 || c6%m6!=0) break; e=ee; }
    ls+=e*log((double)l); }
  return ls; }
static void eval_set(const int *a){
  // p = prod (x-a_i), coefficients p[0..10]
  i128 p[11]={0}; p[0]=1; // p[k] = coefficient of x^k ; build
  for(int i=0;i<11;i++) p[i]=0; p[0]=1; int deg=0;
  for(int i=0;i<10;i++){ for(int k=deg+1;k>=1;k--) p[k]=p[k-1]-a[i]*p[k]; p[0]=-a[i]*p[0]; deg++; }
  // q monic quintic: q5=1
  i128 t;
  if(p[9]&1) return; i128 q4=p[9]/2;
  t=p[8]-q4*q4; if(t&1) return; i128 q3=t/2;
  t=p[7]-2*q4*q3; if(t&1) return; i128 q2=t/2;
  t=p[6]-2*q4*q2-q3*q3; if(t&1) return; i128 q1=t/2;
  t=p[5]-2*q4*q1-2*q3*q2; if(t&1) return; i128 q0=t/2;
  i128 r4=2*q4*q0+2*q3*q1+q2*q2-p[4];
  i128 r3=2*q3*q0+2*q2*q1-p[3];
  i128 r2=2*q2*q0+q1*q1-p[2];
  i128 r1=2*q1*q0-p[1];
  i128 r0=q0*q0-p[0];
  if(r4==0) return;
  { const long long MP[3]={1000000007LL,998244353LL,1000000009LL}; int zero=1;
    for(int k=0;k<3&&zero;k++){ long long m=MP[k];
      #define MM(z) ((long long)(((z)%m+m)%m))
      i128 A_=MM(r4),B_=MM(r3),C_=MM(r2),D_=MM(r1),E_=MM(r0);
      i128 I_=((12*A_%m*E_ - 3*B_%m*D_ + C_*C_)%m+m)%m;
      i128 J_=((72*A_%m*C_%m*E_ + 9*B_%m*C_%m*D_ - 27*A_%m*D_%m*D_ - 27*E_%m*B_%m*B_ - 2*C_%m*C_%m*C_)%m+m)%m;
      i128 Dl=((4*I_%m*I_%m*I_ - J_*J_)%m+m)%m; if(Dl!=0) zero=0; }
    if(zero) return; }
  long double A=r4,B=r3,C=r2,D=r1,E=r0;
  long double I=12*A*E-3*B*D+C*C, J=72*A*C*E+9*B*C*D-27*A*D*D-27*E*B*B-2*C*C*C;
  long double c4=1296.0L*I, c6=23328.0L*J;
  long double m=fmaxl(3*logl(fabsl(c4)+1),2*logl(fabsl(c6)+1)); double H=(double)(m/12);
  // minimal model typically smaller by u in {1,2,3,6,...}: allow slack 2.6
  if(H>HMAX+6.0) return;
  H-=lscale(r4,r3,r2,r1,r0);
  if(H>HMAX) return;
  double S=0, S1v=0;
  for(int i=0;i<NPR;i++){ int pp=PR[i]; if(i==NPR1){ if(S<SMIN) return; S1v=S; } long long R4=(long long)(r4%pp),R3=(long long)(r3%pp),R2=(long long)(r2%pp),R1=(long long)(r1%pp),R0=(long long)(r0%pp);
    R4=(R4+pp)%pp;R3=(R3+pp)%pp;R2=(R2+pp)%pp;R1=(R1+pp)%pp;R0=(R0+pp)%pp;
    if(R4==0) continue; // skip bad-ish
    { long long Ip=((12*R4%pp*R0 - 3*R3%pp*R1 + R2*R2)%pp+pp)%pp;
      long long Jp=((72*R4%pp*R2%pp*R0 + 9*R3%pp*R2%pp*R1 - 27*R4%pp*R1%pp*R1 - 27*R0%pp*R3%pp*R3 - 2*R2%pp*R2%pp*R2)%pp+pp)%pp;
      long long Dp=((4*Ip%pp*Ip%pp*Ip - Jp*Jp)%pp+pp)%pp; if(Dp==0) continue; }
    long long s=0; { long long F[5]; for(int k=0;k<5;k++){ long long x=k; F[k]=((((R4*x+R3)%pp*x+R2)%pp*x+R1)%pp*x+R0)%pp; }
      long long d0=F[0], d1=(F[1]-F[0]+pp)%pp, d2=((F[2]-2*F[1]+F[0])%pp+2*pp)%pp, d3=((F[3]-3*F[2]+3*F[1]-F[0])%pp+4*pp)%pp, d4=((F[4]-4*F[3]+6*F[2]-4*F[1]+F[0])%pp+8*pp)%pp;
      const signed char *ch=CHI[i];
      for(long long x=0;x<pp;x++){ s+=ch[d0]; d0+=d1; if(d0>=pp)d0-=pp; d1+=d2; if(d1>=pp)d1-=pp; d2+=d3; if(d2>=pp)d2-=pp; d3+=d4; if(d3>=pp)d3-=pp; } }
    s+=CHI[i][R4];
    // a_p = -s ; weight -a_p log p/p = s*W
    S+=s*W[i];
    if((i==12 && S<SMIN*0.296) || (i==24 && S<SMIN*0.61)) return; // early abort
  }
  if(S<T2) return;
  #pragma omp critical
  { for(int i=0;i<10;i++) printf("%d ",a[i]); printf("%.4f %.3f %.3f\n",H,S1v,S); fflush(stdout); }
}
static unsigned long long rng_state;
static inline unsigned long long xr(unsigned long long *s){ *s^=*s<<13; *s^=*s>>7; *s^=*s<<17; return *s; }
int main(int argc,char**argv){
  M=atoi(argv[1]); G=atoi(argv[2]); HMAX=atof(argv[3]); SMIN=atof(argv[4]);
  long long NS=argc>5?atoll(argv[5]):1000000; unsigned long long seed=argc>6?atoll(argv[6]):1;
  for(int pp=5;pp<=3000;pp++){int ok=1;for(int d=2;d*d<=pp;d++) if(pp%d==0){ok=0;break;} if(!ok) continue; PR[NPR]=pp; CHI[NPR]=malloc(pp); for(int x=0;x<pp;x++) CHI[NPR][x]=-1; CHI[NPR][0]=0; for(long long x=1;x<pp;x++) CHI[NPR][x*x%pp]=1; W[NPR]=log(pp)/pp; NPR++;}
  for(NPR1=0;NPR1<NPR && PR[NPR1]<=199;NPR1++);
  MODE=getenv("MODE")?atoi(getenv("MODE")):0; T2=getenv("T2")?atof(getenv("T2")):-1e9;
  #pragma omp parallel
  { unsigned long long st=seed*0x9E3779B97F4A7C15ULL+omp_get_thread_num()*7919+1;
    #pragma omp for schedule(dynamic,10000)
    for(long long it=0;it<NS;it++){
      int a[10]; int L[5],g[5];
      // sample 4 distinct sorted centers L2..L5 in [2,M-1], gaps in [1,G]
      int ok=1; int SP=M/4;
      if(MODE==2){ // ~80% roots multiples of k in {4,6}, spread R in [100,150]
        int k=(xr(&st)&1)?4:6; int R=100+xr(&st)%51; int c[9]; int nm=0;
        for(int i=0;i<9;i++){ if((xr(&st)%10)<8){ c[i]=k*(1+xr(&st)%(R/k)); nm++; } else c[i]=1+xr(&st)%R; }
        for(int i=0;i<9;i++) for(int j=i+1;j<9;j++) if(c[j]<c[i]){int t=c[i];c[i]=c[j];c[j]=t;}
        for(int i=0;i<8;i++) if(c[i]==c[i+1]) ok=0; if(!ok) continue; a[0]=0; for(int i=0;i<9;i++) a[i+1]=c[i]; eval_set(a); continue; }
      if(MODE==1){ int c[9]; for(int i=0;i<9;i++) c[i]=1+xr(&st)%M; for(int i=0;i<9;i++) for(int j=i+1;j<9;j++) if(c[j]<c[i]){int t=c[i];c[i]=c[j];c[j]=t;}
        for(int i=0;i<8;i++) if(c[i]==c[i+1]) ok=0; if(!ok) continue; a[0]=0; for(int i=0;i<9;i++) a[i+1]=c[i]; eval_set(a); continue; }
      for(int i=0;i<5;i++) g[i]=1+xr(&st)%G;
      L[0]=0; for(int i=1;i<5;i++){ int sp=1+xr(&st)%SP; L[i]=L[i-1]+g[i-1]+sp; }
      if(L[4]+g[4]>M) ok=0;
      if(!ok) continue;
      for(int i=0;i<5;i++){ a[2*i]=L[i]; a[2*i+1]=L[i]+g[i]; }
      eval_set(a);
    }
  }
  return 0;
}
