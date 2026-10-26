// Z/3 sieve: curves y^2 + a x y + b y = x^3, model [a,0,b,0,0], a>=1.
// usage: z3s Hmax margin thr1 thr2 P1 P2 alo ahi ; env P3 THR3 HLO SLACK
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
typedef long long ll;
static ll pw(ll a, ll e, ll p){ll r=1;a%=p;if(a<0)a+=p;while(e){if(e&1)r=r*a%p;a=a*a%p;e>>=1;}return r;}
static ll inv(ll a, ll p){a%=p;if(a<0)a+=p;if(!a)return -1;return pw(a,p-2,p);}
static int8_t *chi_tab;
static int apf(ll p, ll b2, ll b4, ll b6){ long s=0; for(ll x=0;x<p;x++){ ll f=(4*x%p*x%p*x + b2*x%p*x + 2*b4%p*x + b6)%p; s+=chi_tab[f]; } return (int)(-s); }
static double H3(double a,double b){ double c4=a*a*a*a-24*a*b, c6=-a*a*a*a*a*a+36*a*a*a*b-216*b*b;
  return fmax(3*log(fabs(c4)+1e-300),2*log(fabs(c6)+1e-300))/12; }
static ll gcdl(ll a, ll b){ a=llabs(a); b=llabs(b); while(b){ll r=a%b;a=b;b=r;} return a; }
static int nonmin(ll a, ll b){ ll g=gcdl(a,b); if(g==1) return 0; for(ll l=2;l*l<=g;l++){ if(g%l==0){ if(b%(l*l*l)==0) return 1; while(g%l==0) g/=l; } } if(g>1 && b%(g*g*g)==0) return 1; return 0; }
int main(int argc,char**argv){
  if(argc<9){fprintf(stderr,"usage: z3s Hmax margin thr1 thr2 P1 P2 alo ahi\n");return 1;}
  double Hmax=atof(argv[1]),margin=atof(argv[2]),thr1=atof(argv[3]),thr2=atof(argv[4]);
  int P1=atoi(argv[5]),P2=atoi(argv[6]); ll alo=atoll(argv[7]),ahi=atoll(argv[8]);
  double Hlo=getenv("HLO")?atof(getenv("HLO")):-1e9;
  int P3=getenv("P3")?atoi(getenv("P3")):P2; double thr3=getenv("THR3")?atof(getenv("THR3")):-1e9;
  int maxp=P3; char*sv=calloc(maxp+1,1); int *plist=malloc(sizeof(int)*maxp); int np=0;
  for(int i=2;i<=maxp;i++){ if(!sv[i]){ if(i>=5) plist[np++]=i; for(ll j=(ll)i*i;j<=maxp;j+=i) sv[j]=1;}}
  int n1=0; while(n1<np && plist[n1]<=P1) n1++;
  int n2=n1; while(n2<np && plist[n2]<=P2) n2++;
  float **W=malloc(sizeof(float*)*np), **W0=malloc(sizeof(float*)*np);
  char tfn[128]; sprintf(tfn,"tab3_%d.bin",P3); FILE*tf=fopen(tfn,"rb"); int loaded=0;
  if(tf){ loaded=1; for(int i=0;i<np;i++){ int p=plist[i]; W[i]=malloc(sizeof(float)*p); W0[i]=malloc(sizeof(float)*p);
      if(fread(W[i],sizeof(float),p,tf)!=(size_t)p || fread(W0[i],sizeof(float),p,tf)!=(size_t)p) loaded=0; } fclose(tf); }
  if(!loaded) for(int i=0;i<np;i++){ ll p=plist[i];
    chi_tab=malloc(p); memset(chi_tab,-1,p); chi_tab[0]=0; for(ll x=1;x<p;x++) chi_tab[x*x%p]=1;
    W[i]=malloc(sizeof(float)*p); W0[i]=malloc(sizeof(float)*p); double lw=log((double)p)/p;
    ll i27=inv(27,p);
    for(ll s=0;s<p;s++){ // curve [1,0,s,0,0]: b2=1,b4=s,b6=s^2 ; singular if s=0 or s=1/27
      if(s==0 || s==i27){W[i][s]=0;} else W[i][s]=(float)(-apf(p,1,s,s*s%p)*lw);
      // a==0 mod p: curve [0,0,j,0,0]: b2=0,b4=0,b6=j^2
      if(s==0) W0[i][s]=0; else W0[i][s]=(float)(-apf(p,0,0,s*s%p)*lw); }
    free(chi_tab);}
  if(!loaded && getenv("SAVETAB")){ tf=fopen(tfn,"wb"); for(int i=0;i<np;i++){ fwrite(W[i],sizeof(float),plist[i],tf); fwrite(W0[i],sizeof(float),plist[i],tf);} fclose(tf); }
  const int BL=2048;
  int8_t *S8=malloc(BL); int8_t **Q8=malloc(sizeof(int8_t*)*n1); for(int i=0;i<n1;i++) Q8[i]=malloc(plist[i]+BL);
  double slack=getenv("SLACK")?atof(getenv("SLACK")):2.0; int8_t ithr8=(int8_t)lrint(fmax(-127,fmin(127,(thr1-slack)*4)));
  int bmode=getenv("BMODE")!=NULL;
  double Hcut=Hmax+margin; ll cnt=0,pass1=0,tot=0;
  ll *ia3=malloc(sizeof(ll)*np);
  for(ll a=alo; a<ahi; a++){
    double da=(double)a;
    if(log(da)>Hcut+0.05) continue;
    double B=1; while(B<4e18 && (H3(da,B)<Hcut+0.2 || H3(da,-B)<Hcut+0.2 || B<da*da*da)) B*=1.25;
    ll blo=-(ll)B, bhi=(ll)B;
    if(getenv("SLO")){ double lo=atof(getenv("SLO"))*da*da*da; if(lo>-B) blo=(ll)floor(lo); }
    if(getenv("SHI")){ double hi=atof(getenv("SHI"))*da*da*da; if(hi<B) bhi=(ll)ceil(hi); }
    if(blo>bhi) continue;
    ll BA=1, BC=0;
    if(bmode){ BA=6; ll c3=((-a)%3+3)%3, c2=((a%2)+2)%2; for(ll c=0;c<6;c++) if(c%3==c3 && c%2==c2){ BC=c; break; } }
    for(int i=0;i<np;i++){ ll p=plist[i]; ll am=((a%p)+p)%p; ia3[i]= am? inv(am*am%p*am%p,p) : -1; }
    for(int i=0;i<n1;i++){ ll p=plist[i]; int8_t*q=Q8[i];
      { ll bm=BC%p, st=BA%p; for(int j=0;j<p;j++){ q[j]=(int8_t)lrint(((ia3[i]<0)? W0[i][bm] : W[i][bm*ia3[i]%p])*4); bm+=st; if(bm>=p) bm-=p; } }
      for(int j=p;j<p+BL;j++) q[j]=q[j-p]; }
    ll wlo=(ll)floor((double)(blo-BC)/BA), whi=(ll)floor((double)(bhi-BC)/BA);
    for(ll w0=wlo; w0<=whi; w0+=BL){
      int len=(whi-w0+1<BL)?(int)(whi-w0+1):BL; ll b0=BA*w0+BC;
      double h0=H3(da,(double)b0), h1=H3(da,(double)b0+BA*(len/2)), h2=H3(da,(double)b0+BA*(len-1));
      if(fmin(h0,fmin(h1,h2))>Hcut+0.02) continue;
      if(fmax(h0,fmax(h1,h2))<Hlo-0.02) continue;
      tot+=len;
      memset(S8,0,BL);
      for(int i=0;i<n1;i++){ int p=plist[i]; ll off=((w0%p)+p)%p; const int8_t*q=Q8[i]+off; for(int j=0;j<BL;j++) S8[j]+=q[j]; }
      for(int j=0;j<len;j++){ if(S8[j]<ithr8) continue; ll b=b0+BA*j; if(b==0) continue;
        double s1=0; for(int i=0;i<n1;i++){ ll p=plist[i]; ll bm=((b%p)+p)%p; s1+= (ia3[i]<0)? W0[i][bm] : W[i][bm*ia3[i]%p]; }
        if(s1<thr1) continue;
        if(nonmin(a,b)) continue;
        pass1++;
        double h=H3(da,(double)b); if(h>Hcut || h<Hlo) continue;
        double s2=s1;
        for(int i=n1;i<n2;i++){ ll p=plist[i]; ll bm=((b%p)+p)%p; s2+= (ia3[i]<0)? W0[i][bm] : W[i][bm*ia3[i]%p]; }
        if(s2<thr2) continue;
        double s3=s2;
        for(int i=n2;i<np;i++){ ll p=plist[i]; ll bm=((b%p)+p)%p; s3+= (ia3[i]<0)? W0[i][bm] : W[i][bm*ia3[i]%p]; }
        if(s3<thr3) continue;
        printf("%lld %lld %.3f %.3f %.4f %.3f\n",a,b,s1,s2,h,s3); cnt++; }
    }
  }
  fflush(stdout);
  fprintf(stderr,"done a=[%lld,%lld) tot=%lld pass1=%lld out=%lld\n",alo,ahi,tot,pass1,cnt);
  return 0;
}
