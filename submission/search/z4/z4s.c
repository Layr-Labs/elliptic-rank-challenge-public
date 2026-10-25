// One-parameter torsion family sieve: E_t Tate normal form, t=u/v.
// usage: tsieve N T Hmax margin thr1 thr2 P1 P2 vlo vhi [ulo uhi]
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
typedef long long ll;
static int N;
static ll pw(ll a, ll e, ll p){ll r=1;a%=p;if(a<0)a+=p;while(e){if(e&1)r=r*a%p;a=a*a%p;e>>=1;}return r;}
static ll inv(ll a, ll p){a%=p;if(a<0)a+=p;if(!a)return -1;return pw(a,p-2,p);}
// compute b,c mod p from t mod p; return 0 on pole
static int fam_modp(ll t, ll p, ll *b, ll *c){
  ll d,m,f,i;
  switch(N){
  case 4: *b=t;*c=0;return 1;
  case 5: *b=t;*c=t;return 1;
  case 6: *b=(t+t*t)%p;*c=t;return 1;
  case 7: *b=((t*t%p)*t - t*t%p + p)%p; *c=(t*t-t+p)%p; return 1;
  case 8: *b=(2*t-1+p)%p*((t-1+p)%p)%p; i=inv(t,p); if(i<0)return 0; *c=*b*i%p; return 1;
  case 9: *c=t*t%p*((t-1+p)%p)%p; *b=*c*((t*t-t+1+p)%p)%p; return 1;
  case 10:{ ll den=(t-(t-1)*(t-1))%p; den=(den%p+p)%p; i=inv(den,p); if(i<0)return 0; d=t*t%p*i%p; *c=((t*d-t)%p+p)%p; *b=*c*d%p; return 1;}
  case 12:{ i=inv(t-1,p); if(i<0) return 0; m=((3*t-3*t%p*t-1)%p+p)%p*i%p; ll j=inv(1-t,p); f=m*j%p; d=(m+t)%p; *c=f*((d-1+p)%p)%p; *b=*c*d%p; return 1;}
  case 24:{ ll i16=inv(16,p); *b=((t*t-i16)%p+p)%p; *c=0; return 1;}
  case 26:{ i=inv(t*t-9,p); if(i<0) return 0; *c=((10-2*t)%p+p)%p*i%p; *b=(*c+*c**c)%p; return 1;}
  case 28:{ ll U=t*((8*t+2)%p)%p, V=((8*t%p*t-1)%p+p)%p; i=inv(V,p); if(i<0) return 0; ll tt=U*i%p; int sv=N; N=8; int r=fam_modp(tt,p,b,c); N=sv; return r;}
  }
  return 0;
}
static int fam_d(double t, double *b, double *c){
  double d,m,f;
  switch(N){
  case 4: *b=t;*c=0;return 1;
  case 5: *b=t;*c=t;return 1;
  case 6: *b=t+t*t;*c=t;return 1;
  case 7: *b=t*t*t-t*t;*c=t*t-t;return 1;
  case 8: *b=(2*t-1)*(t-1); if(t==0)return 0; *c=*b/t; return 1;
  case 9: *c=t*t*(t-1); *b=*c*(t*t-t+1); return 1;
  case 10: d=t*t/(t-(t-1)*(t-1)); *c=t*d-t; *b=*c*d; return 1;
  case 12: m=(3*t-3*t*t-1)/(t-1); f=m/(1-t); d=m+t; *c=f*(d-1); *b=*c*d; return 1;
  case 24: *b=t*t-1.0/16; *c=0; return 1;
  case 26: *c=(10-2*t)/(t*t-9); *b=*c+*c**c; return 1;
  }
  return 0;
}
static double lam_d(double u,double v){
  switch(N){case 4:case 5:case 6:return v;case 7:return v*v;case 8:return u*v;case 9:return v*v*v;
  case 10:return v*(u*u-3*u*v+v*v);case 12:return v*(u-v)*(u-v)*(u-v);case 24:return 4*v;case 26:return (u-3*v)*(u+3*v);}
  return 1;
}
static double mlog[64][64]; static int have_m=0;
static void load_m(void){ char fn[64]; sprintf(fn,"mtab%d.plain",N); FILE*f=fopen(fn,"r"); if(!f) return;
  for(int i=0;i<64;i++) for(int j=0;j<64;j++){ int m; if(fscanf(f,"%d",&m)!=1){fclose(f);return;} mlog[i][j]=(m>0)?log((double)m):0; }
  fclose(f); have_m=1; }
static double Hint0(ll u, ll v);
// returns log-height H of model after the known 2-adic reduction
static double Hint(ll u, ll v){ double h=Hint0(u,v); if(have_m){ h-=mlog[((u%64)+64)%64][((v%64)+64)%64]; } return h; }
static double Hint0(ll u, ll v){
  if(N==28){ ll U=u*(8*u+2*v), V=8*u*u-v*v; ll a=llabs(U),bb=llabs(V); while(bb){ll r=a%bb;a=bb;bb=r;} if(a==0) return 1e9; U/=a; V/=a; if(V<0){V=-V;U=-U;} N=8; double h=Hint(U,V); N=28; return h; }
  double t=(double)u/(double)v,b,c; if(!fam_d(t,&b,&c)) return 1e9;
  double L=fabs(lam_d(u,v)); if(L==0) return 1e9;
  double a1=1-c,a2=-b,a3=-b;
  double b2=a1*a1+4*a2,b4=a1*a3,b6=a3*a3;
  double c4=b2*b2-24*b4, c6=-b2*b2*b2+36*b2*b4-216*b6;
  double x=fmax(3*log(fabs(c4)+1e-300), 2*log(fabs(c6)+1e-300));
  return x/12+log(L);
}
static int *plist; static int np;
static int8_t *chi_tab;
static int ap_count(ll p, ll b, ll c){ // a_p for y^2+(1-c)xy-by=x^3-bx^2
  ll a1=(1-c+p)%p,a2=(p-b)%p,a3=(p-b)%p;
  ll b2=(a1*a1+4*a2)%p,b4=a1*a3%p,b6=a3*a3%p;
  // sum chi(4x^3+b2x^2+2b4x+b6)
  long s=0;
  for(ll x=0;x<p;x++){ ll f=((4*x%p*x%p*x + b2*x%p*x + 2*b4*x + b6)%p); s+=chi_tab[f]; }
  return (int)(-s);
}

// Z/4 sieve: t=u/(D s^2), model [D s, -D u, -D^2 s u, 0, 0]
static double H4(double D,double s,double u){
  double b2=D*D*s*s-4*D*u, b4=-D*D*D*s*s*u, b6=D*D*D*D*s*s*u*u;
  double c4=b2*b2-24*b4, c6=-b2*b2*b2+36*b2*b4-216*b6;
  return fmax(3*log(fabs(c4)+1e-300),2*log(fabs(c6)+1e-300))/12; }
int main(int argc,char**argv){
  if(argc<10){fprintf(stderr,"usage: z4s Hmax margin thr1 thr2 P1 P2 D slo shi [step off]\n");return 1;}
  N=4; double Hmax=atof(argv[1]),margin=atof(argv[2]),thr1=atof(argv[3]),thr2=atof(argv[4]);
  int P1=atoi(argv[5]),P2=atoi(argv[6]); ll D=atoll(argv[7]); ll slo=atoll(argv[8]),shi=atoll(argv[9]);
  ll step=1,off=0; if(argc>11){step=atoll(argv[10]);off=atoll(argv[11]);}
  double Hlo=getenv("HLO")?atof(getenv("HLO")):-1e9;
  int P3=getenv("P3")?atoi(getenv("P3")):P2; double thr3=getenv("THR3")?atof(getenv("THR3")):-1e9;
  int maxp=P3; char*sv=calloc(maxp+1,1); plist=malloc(sizeof(int)*maxp); np=0;
  for(int i=2;i<=maxp;i++){ if(!sv[i]){ if(i>=3) plist[np++]=i; for(ll j=(ll)i*i;j<=maxp;j+=i) sv[j]=1;}}
  int n1=0; while(n1<np && plist[n1]<=P1) n1++;
  int n2=n1; while(n2<np && plist[n2]<=P2) n2++;
  int P0=getenv("P0")?atoi(getenv("P0")):P1; int n0=0; while(n0<n1 && plist[n0]<=P0) n0++;
  short ithr0=getenv("THR0")?(short)lrint(atof(getenv("THR0"))*64):-32000;
  float **W=malloc(sizeof(float*)*np);
  char tfn[128]; sprintf(tfn,"tab%d_%d.bin",N,P3); FILE*tf=fopen(tfn,"rb"); int loaded=0;
  if(tf){ loaded=1; for(int i=0;i<np;i++){ int p=plist[i]; W[i]=malloc(sizeof(float)*(p+1)); if(fread(W[i],sizeof(float),p+1,tf)!=(size_t)(p+1)) loaded=0; } fclose(tf); }
  if(!loaded) for(int i=0;i<np;i++){ int p=plist[i];
    chi_tab=malloc(p); memset(chi_tab,-1,p); chi_tab[0]=0; for(ll x=1;x<p;x++) chi_tab[x*x%p]=1;
    W[i]=malloc(sizeof(float)*(p+1)); double lw=log((double)p)/p;
    for(ll t=0;t<p;t++){ ll b,c; if(!fam_modp(t,p,&b,&c)) {W[i][t]=0;continue;} int a=ap_count(p,b,c); W[i][t]=(float)(-a*lw);}
    W[i][p]=0; free(chi_tab);}
  if(!loaded && getenv("SAVETAB")){ tf=fopen(tfn,"wb"); for(int i=0;i<np;i++) fwrite(W[i],sizeof(float),plist[i]+1,tf); fclose(tf); }
  const int BL=2048; short *S=malloc(sizeof(short)*BL); short **Q=malloc(sizeof(short*)*n1); short ithr1=(short)lrint(thr1*64);
  int8_t *S8=malloc(BL); int8_t **Q8=malloc(sizeof(int8_t*)*n1); for(int i=0;i<n1;i++) Q8[i]=malloc(plist[i]+BL);
  double slack=getenv("SLACK")?atof(getenv("SLACK")):2.0; int8_t ithr8=(int8_t)lrint(fmax(-127,fmin(127,(thr1-slack)*4)));
  ll UA=getenv("UA")?atoll(getenv("UA")):1, UB=getenv("UB")?atoll(getenv("UB")):0;
  int skip2=getenv("SKIP2")!=NULL;
  double Hcut=Hmax+margin; ll cnt=0,pass1=0,tot=0;
  ll *ivs=malloc(sizeof(ll)*np);
  for(ll s=slo; s<shi; s++){ if(((s-off)%step+step)%step) continue; if(skip2 && s%4==2) continue;
    double ds=(double)s, dD=(double)D;
    { double mh=H4(dD,ds,0.5); for(int k=-600;k<=600;k++){ double tt=k/200.0; double h=H4(dD,ds,tt*dD*ds*ds); if(h<mh) mh=h; } if(mh>Hcut+0.3) continue; }
    // find u range
    double U=1; while(U<4e18 && (H4(dD,ds,U)<Hcut+0.2 || H4(dD,ds,-U)<Hcut+0.2 || U<4*dD*ds*ds)) U*=1.25;
    if(U>4e18) U=4e18;
    ll ulo=-(ll)U, uhi=(ll)U;
    if(getenv("TLO")){ double lo=atof(getenv("TLO"))*dD*ds*ds; if(lo>-U) ulo=(ll)floor(lo); }
    if(getenv("THI")){ double hi=atof(getenv("THI"))*dD*ds*ds; if(hi<U) uhi=(ll)ceil(hi); }
    if(ulo>uhi) continue;
    ll v=D*s*s; // v up to ~1e18 OK for D*s^2 < 9e18
    for(int i=0;i<np;i++) ivs[i]=inv(v%plist[i],plist[i]);
    for(int i=0;i<n1;i++){ int p=plist[i]; ll iv=inv(v%p,p); int8_t*q=Q8[i];
      if(iv<0){ for(int j=0;j<p+BL;j++) q[j]=(int8_t)lrint(W[i][p]*4); }
      else { ll st=UA%p*iv%p, t=UB%p*iv%p; for(int j=0;j<p;j++){ q[j]=(int8_t)lrint(W[i][t]*4); t+=st; if(t>=p)t-=p; } for(int j=p;j<p+BL;j++) q[j]=q[j-p]; } }
    ll wlo=(ll)floor((double)(ulo-UB)/UA), whi=(ll)floor((double)(uhi-UB)/UA);
    for(ll w0=wlo; w0<=whi; w0+=BL){
      int len = (whi-w0+1<BL)?(int)(whi-w0+1):BL; ll u0=UA*w0+UB;
      double h0=H4(dD,ds,(double)u0), h1=H4(dD,ds,(double)u0+UA*(len/2)), h2=H4(dD,ds,(double)u0+UA*(len-1));
      double mh=fmin(h0,fmin(h1,h2)); if(mh>Hcut+0.02) continue;
      double xh=fmax(h0,fmax(h1,h2)); if(xh<Hlo-0.02) continue;
      tot+=len;
            memset(S8,0,BL);
      for(int i=0;i<n1;i++){ int p=plist[i]; ll off2=((w0%p)+p)%p; const int8_t*q=Q8[i]+off2; for(int j=0;j<BL;j++) S8[j]+=q[j]; }
      for(int j=0;j<len;j++){ if(S8[j]<ithr8) continue; ll u=u0+UA*j; if(u==0) continue;
        double s1=0; for(int i=0;i<n1;i++){ int p=plist[i]; ll iv=ivs[i]; ll t; if(iv<0) t=p; else { t=((u%p)+p)%p*iv%p; } s1+=W[i][t]; }
        if(s1<thr1) continue; S[j]=(short)lrint(s1*64);
        ll a=llabs(u),bb=D*s; while(bb){ll r=a%bb;a=bb;bb=r;} if(a!=1) continue;
        pass1++;
        double h=H4(dD,ds,(double)u); if(h>Hcut || h<Hlo) continue;
        double s2=S[j]/64.0;
        for(int i=n1;i<n2;i++){ int p=plist[i]; ll iv=ivs[i]; ll t; if(iv<0) t=p; else { t=((u%p)+p)%p*iv%p; } s2+=W[i][t]; }
        if(s2<thr2) continue;
        double s3=s2;
        for(int i=n2;i<np;i++){ int p=plist[i]; ll iv=ivs[i]; ll t; if(iv<0) t=p; else { t=((u%p)+p)%p*iv%p; } s3+=W[i][t]; }
        if(s3<thr3) continue;
        if(np>n2) printf("%lld %lld %lld %.3f %.3f %.4f %.3f\n",u,D,s,S[j]/64.0,s2,h,s3); else printf("%lld %lld %lld %.3f %.3f %.4f\n",u,D,s,S[j]/64.0,s2,h); cnt++; }
    }
  }
  fflush(stdout);
  fprintf(stderr,"done D=%lld s=[%lld,%lld) tot=%lld pass1=%lld out=%lld\n",D,slo,shi,tot,pass1,cnt);
  return 0;
}
