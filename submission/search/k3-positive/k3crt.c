// Original CRT enumeration and scaled-height shortlist for Elkies Z/4 K3.
// The finite-field table routines are reused from zereraz's inherited kub.c.
#define KUB_NO_MAIN
#include "../gpu-table-sieve/src/kub.c"
#include <omp.h>
typedef __int128 wide;
typedef struct {long long mod,res,oddscale;int vtwo;long long cap;} workcase;
static workcase jobs[4000];static int njob;
static unsigned long long positions=0,eligible=0,hits=0;
static int vp(wide x,int p){if(x<0)x=-x;if(!x)return 99;int e=0;while(x%p==0){x/=p;e++;}return e;}
static wide aval(long long p,long long q){return ((wide)8*p-q)*((wide)32*p+7*q);}
static wide bval(long long p,long long q){return (wide)8*q*(p+q)*((wide)15*p-8*q)*((wide)31*p-7*q);}
static double scaledheight(long long p,long long q){wide aa=aval(p,q),bb=bval(p,q);if(!aa||!bb||aa*aa==16*bb)return 1e100;long double u=1;int ps[]={2,3,5,7,19};for(int k=0;k<5;k++){int l=ps[k],va=vp(aa,l),vb=vp(bb,l)/2,e=va<vb?va:vb;while(e--)u*=l;}long double a=(long double)aa/u,b=(long double)bb/(u*u),c4=a*a*a*a-16*a*a*b+16*b*b,c6=-a*a*a*a*a*a+24*a*a*a*a*b-120*a*a*b*b-64*b*b*b;return fmaxl(logl(fabsl(c4))/4,logl(fabsl(c6))/6)-log(2.0);}
static long long invll(long long a,long long m){long long t=0,nt=1,r=m,nr=a%m;while(nr){long long q=r/nr,z=t-q*nt;t=nt;nt=z;z=r-q*nr;r=nr;nr=z;}return (t%m+m)%m;}
static void addclasses(int index,long long mod,long long residue,long long scale,double H){int lp[]={3,5,7,19};if(index==4){for(int kt=0;kt<=13;kt++){long long scale2=kt==0?1:1LL<<(((3+kt)/2)<8?((3+kt)/2):8);if(kt==1)scale2=4;long double bound=expl((H+logl((long double)scale*scale2*2))/2);if(bound>30000000)bound=30000000;jobs[njob++]=(workcase){mod,residue,scale,kt,(long long)ceill(bound)};if(njob>=4000){fprintf(stderr,"case capacity\n");exit(2);}}return;}
 int l=lp[index],maxe=l==19?1:2;addclasses(index+1,mod,residue,scale,H);
 for(int e=1;e<=maxe;e++){long long pe=1,me=1;for(int i=0;i<e;i++)pe*=l;me=pe*pe;for(long long r=0;r<me;r++){if(vp(aval(r,1),l)<e||vp(bval(r,1),l)<2*e)continue;long long crt=(residue+mod*(((r-residue)%me+me)%me*invll(mod,me)%me))%(mod*me);addclasses(index+1,mod*me,crt,scale*pe,H);}}
}
int main(int argc,char**argv){if(argc<3){fprintf(stderr,"usage: k3crt Hmax threads [S300 S2000 S10000]\n");return 1;}double hmax=atof(argv[1]);int nth=atoi(argv[2]);double t1=argc>3?atof(argv[3]):10,t2=argc>4?atof(argv[4]):24,t3=argc>5?atof(argv[5]):38;FAM=404;load_family(404);build_tables(nth);addclasses(0,1,0,1,hmax);fprintf(stderr,"cases=%d H<=%.6f thresholds=%.2f/%.2f/%.2f cap=30000000\n",njob,hmax,t1,t2,t3);double began=now();
 #pragma omp parallel for schedule(dynamic,1) num_threads(nth)
 for(int job=0;job<njob;job++){workcase w=jobs[job];unsigned long long n=0,ne=0,nh=0;long long qstep=w.vtwo==13?8192:(1LL<<(w.vtwo+1)),qstart=w.vtwo==0?1:1LL<<w.vtwo;
 for(long long q=qstart;q<=w.cap;q+=qstep){if(gcd64(q,w.mod)!=1)continue;long long rr=(w.res*(q%w.mod))%w.mod;long long p0=-w.cap+(((rr+w.cap)%w.mod+w.mod)%w.mod);for(long long p=p0;p<=w.cap;p+=w.mod){n++;if(!p||gcd64(p,q)!=1)continue;if(w.vtwo&&!(p&1))continue;double H=scaledheight(p,q);if(H>=hmax)continue;ne++;double s1=0,s2=0,s3=0;int sv=0;for(int i=0;i<n1;i++){long long l=primes[i],qm=q%l;if(qm){long long pm=(p%l+l)%l;sv+=tab[i][pm*inv[i][qm]%l];}}if(sv/SCALE<t1)continue;if(!sums_staged(p,q,t2,t3,&s1,&s2,&s3))continue;if(s3<t3||s1<t1)continue;nh++;
 #pragma omp critical
 {printf("%lld %lld %.2f %.2f %.2f %.6f\n",p,q,s1,s2,s3,H);fflush(stdout);}
 }}
 #pragma omp atomic
 positions+=n;
 #pragma omp atomic
 eligible+=ne;
 #pragma omp atomic
 hits+=nh;
 }
 fprintf(stderr,"done positions=%llu height_eligible=%llu output=%llu seconds=%.3f (domains overlap; height allows extra 2-scale)\n",positions,eligible,hits,now()-began);return 0;}
