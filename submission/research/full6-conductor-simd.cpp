#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <queue>
#include <sstream>
#include <string>
#include <vector>
#include <tuple>
#include <omp.h>
constexpr int SCALE=128, BLOCK=32768, KEEP_LOCAL=20000, KEEP_FINAL=12000;
struct Table {int p;std::vector<int16_t>w;std::vector<int>inv;};
struct Candidate{int p,q,r,score,key;double h,logN;bool operator<(const Candidate&b)const{return key<b.key;}};
using Heap=std::priority_queue<Candidate>;
int main(int argc,char**argv){
 int N=argc>1?std::stoi(argv[1]):150000,threads=argc>2?std::stoi(argv[2]):8;
 std::ifstream data("/tmp/ec-rank-full6-prime-tables.txt");std::string line;std::vector<Table>tables;
 int lower=0,upper=0;
 while(std::getline(data,line)){
  for(char&c:line)if(c=='['||c==']'||c==',')c=' ';
  std::istringstream row(line);Table t;row>>t.p;int a;
  while(row>>a)t.w.push_back(std::lround(SCALE*(a==9999?0.:a*std::log(t.p)/t.p)));
  if(t.w.size()!=t.p)return 2;
  lower+=std::min<int>(0,*std::min_element(t.w.begin(),t.w.end()));
  upper+=std::max<int>(0,*std::max_element(t.w.begin(),t.w.end()));
  t.inv.resize(t.p);t.inv[1]=1;for(int i=2;i<t.p;i++)t.inv[i]=t.p-(t.p/i)*t.inv[t.p%i]%t.p;
  tables.push_back(std::move(t));
 }
 if(lower<-32768||upper>32767){std::cerr<<"COUNTER_RANGE_UNSAFE "<<lower<<" "<<upper<<"\n";return 3;}
 std::cerr<<"PRIMES "<<tables.size()<<" COUNTER_BOUNDS "<<lower<<" "<<upper<<"\n";
 std::vector<double>radlog(8*N+1,0.);std::vector<unsigned char>smallbad(8*N+1,0);
 std::vector<int>marked(8*N+1,0);
 int localU[128][128]={{}},localF[128][128]={{}};
 {std::ifstream local("/tmp/ec-rank-full6-local-tables.txt");int r,s,u,f;while(local>>r>>s>>u>>f){localU[r][s]=u;localF[r][s]=f;}}

 for(int prime=2;prime<=8*N;prime++)if(!marked[prime]){
  double lp=std::log(prime);
  for(int n=prime;n<=8*N;n+=prime){marked[n]=prime;if(prime>3){radlog[n]+=lp;if(prime<=2000)smallbad[n]++;}}
 }
 std::vector<std::vector<Heap>>heaps(threads,std::vector<Heap>(16));
 double bases[]={0,3.560391,5.796933,8.235200,11.329331,14.780907,22.598355};
 omp_set_num_threads(threads);double started=omp_get_wtime();long long counterPairs=0;
 #pragma omp parallel reduction(+:counterPairs)
 {
  int tid=omp_get_thread_num();std::vector<int16_t>counts(BLOCK),period(2001);
  #pragma omp for schedule(dynamic,16)
  for(int q=1;q<=N;q++){
   for(int lo=-N;lo<=N;lo+=BLOCK){
    int len=std::min(BLOCK,N-lo+1);counterPairs+=len;std::fill(counts.begin(),counts.begin()+len,0);
    for(const Table&t:tables){
     int prime=t.p,qr=q%prime;
     if(!qr)continue;
     int inv=t.inv[qr],residue=(int)(((1LL*lo*inv)%prime+prime)%prime);
     for(int i=0;i<prime;i++){period[i]=t.w[residue];residue+=inv;if(residue>=prime)residue-=prime;}
     int start=0;
     for(;start+prime<=len;start+=prime){
      #pragma omp simd
      for(int i=0;i<prime;i++)counts[start+i]+=period[i];
     }
     #pragma omp simd
     for(int i=0;i<len-start;i++)counts[start+i]+=period[i];
    }
    for(int j=0;j<len;j++){
     int score=counts[j];if(score>-10*SCALE)continue;
     int p=lo+j;if(std::gcd(p,q)!=1)continue;
     long long l0=q,l1=1LL*p-q,l2=1LL*p+q,l3=3LL*p-5LL*q,l4=1LL*p-3LL*q,l5=3LL*p-q;
     if(!l1||!l2||!l3||!l4||!l5)continue;
     auto current=std::make_tuple(std::max(std::llabs((long long)p),(long long)q),(long long)q,(long long)p);
     bool representative=true;
     long long images[5][2]={{2LL*q-p,q},{7LL*q-p,5LL*q-3LL*p},{7LL*q-5LL*p,q-3LL*p},{3LL*q-5LL*p,5LL*q-3LL*p},{-5LL*q-p,q-3LL*p}};
     for(auto&image:images){
      long long num=image[0],den=image[1];if(!den)continue;
      long long g=std::gcd(num,den);num/=g;den/=g;if(den<0){num=-num;den=-den;}
      if(std::make_tuple(std::max(std::llabs(num),den),den,num)<current){representative=false;break;}
     }
     if(!representative)continue;

     int rp=(p%128+128)%128,rq=q%128,u2=localU[rp][rq],f2=localF[rp][rq],f3=q%3==0?2:1;
     long double P=p,Q=q,P2=P*P,Q2=Q*Q;
     long double A=37*Q2*Q2-84*P*Q*Q2+102*P2*Q2-36*P2*P*Q-3*P2*P2;
     long double B=32*Q*l1*l1*l1*l2*l2*l2*l3;
     long double c4=16*(A*A-3*B),c6=-64*A*A*A+288*A*B;
     double h=std::max(logl(fabsl(c4))/4,logl(fabsl(c6))/6)-u2*std::log(2.);
     double hr=std::round(h*1e6)/1e6;
     int r=1;while(r<6&&hr>=bases[r])r++;
     if(hr>=bases[6]||r<4)continue;
     double lnN=f2*std::log(2.)+f3*std::log(3.);int nb=2;
     for(long long factor:{l0,l1,l2,l3,l4,l5}){auto v=std::llabs(factor);lnN+=radlog[v];nb+=smallbad[v];}
     for(int mode=0;mode<2;mode++){
      int key=score-std::lround(SCALE*(mode==0?0.4*nb:0.5*lnN));
      auto&heap=heaps[tid][2*r+mode];Candidate c={p,q,r,score,key,h,lnN};
      if(heap.size()<KEEP_LOCAL)heap.push(c);
      else if(key<heap.top().key){heap.pop();heap.push(c);}
     }
    }
   }
   if(q%20000==0){
    #pragma omp critical
    std::cerr<<"DENOM "<<q<<" SECONDS "<<omp_get_wtime()-started<<"\n";
   }
  }
 }
 std::cout<<std::setprecision(17);
 for(int r=4;r<=6;r++)for(int mode=0;mode<2;mode++){
  Heap merged;
  for(auto&byrank:heaps){auto&h=byrank[2*r+mode];while(!h.empty()){auto c=h.top();h.pop();if(merged.size()<KEEP_FINAL)merged.push(c);else if(c.key<merged.top().key){merged.pop();merged.push(c);}}}
  while(!merged.empty()){auto c=merged.top();merged.pop();std::cout<<"62 "<<c.p<<" "<<c.q<<" "<<c.h<<" "<<double(c.score)/SCALE<<" "<<c.logN<<" "<<mode<<"\n";}
 }
 std::cerr<<"FINISHED_GRID "<<counterPairs<<" SECONDS "<<omp_get_wtime()-started<<"\n";
}
