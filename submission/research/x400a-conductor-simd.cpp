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
constexpr int SCALE=128, BLOCK=32768, KEEP_LOCAL=5000, KEEP_FINAL=3000;
struct Table {int p;std::vector<int16_t>w;std::vector<int>inv;};
struct Candidate{int p,q,r,score,key;double h,logN;bool operator<(const Candidate&b)const{return key<b.key;}};
using Heap=std::priority_queue<Candidate>;
int main(int argc,char**argv){
 int N=argc>1?std::stoi(argv[1]):10000,threads=argc>2?std::stoi(argv[2]):8;
 std::ifstream data("/tmp/ec-rank-x400a-prime-tables.txt");std::string line;std::vector<Table>tables;
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

 int bound=225*N+1;
 std::vector<double>radlog(bound,0.);std::vector<unsigned char>smallbad(bound,0);
 std::vector<int>marked(bound,0);
 for(int prime=2;prime<bound;prime++)if(!marked[prime]){
  double lp=std::log(prime);
  for(int n=prime;n<bound;n+=prime){marked[n]=prime;radlog[n]+=lp;if(prime<=2000)smallbad[n]++;}
 }
 std::vector<std::vector<Heap>>heaps(threads,std::vector<Heap>(34));
 double bases[]={0,1.429257,2.033544,2.788376,3.495927,4.628728,5.785321,6.924776,8.108753,9.460889,11.179460,14.601955,18.004452,26.418419,28.082546,33.877449,1e9};
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

     long long ls[12]={q,121LL*q-99LL*p,26LL*p-23LL*q,75LL*p-44LL*q,p,21LL*q-74LL*p,25LL*p-36LL*q,1LL*p+77LL*q,1LL*p-23LL*q,1LL*p-4LL*q,25LL*p-11LL*q,74LL*p-121LL*q};
     bool singular=false;for(auto x:ls)if(!x)singular=true;if(singular)continue;
     long double A=1,B=1;
     for(int k=0;k<4;k++)A*=ls[k];
     for(int k=4;k<8;k++)B*=ls[k];
     long double c4=16*(A*A-A*B+B*B),c6=-64*(A+B)*(A+B)*(A+B)+288*(A+B)*A*B;
     double h=std::max(logl(fabsl(c4))/4,logl(fabsl(c6))/6)-std::log(2.);
     for(int pr:{2,3,5,7,11,17,19,23,41}){
      int va=0,vb=0;
      for(int k=0;k<8;k++){long long v=std::llabs(ls[k]);int val=0;while(v%pr==0){v/=pr;val++;}if(k<4)va+=val;else vb+=val;}
      h-=(std::min(va,vb)/2)*std::log(pr);
     }
     double hr=std::round(h*1e6)/1e6;
     int r=1;while(r<16&&hr>=bases[r])r++;
     if(r<9||r>14)continue;
     double lnN=0;int nb=0;
     for(auto factor:ls){int v=std::llabs(factor);lnN+=radlog[v];nb+=smallbad[v];}
     for(int mode=0;mode<2;mode++){
      int key=score-std::lround(SCALE*(mode==0?0.4*nb:0.5*lnN));
      auto&heap=heaps[tid][2*r+mode];Candidate c={p,q,r,score,key,h,lnN};
      if(heap.size()<KEEP_LOCAL)heap.push(c);
      else if(key<heap.top().key){heap.pop();heap.push(c);}
     }
    }
   }
   if(q%5000==0){
    #pragma omp critical
    std::cerr<<"DENOM "<<q<<" SECONDS "<<omp_get_wtime()-started<<"\n";
   }
  }
 }
 std::cout<<std::setprecision(17);
 for(int r=9;r<=14;r++)for(int mode=0;mode<2;mode++){
  Heap merged;
  for(auto&byrank:heaps){auto&h=byrank[2*r+mode];while(!h.empty()){auto c=h.top();h.pop();if(merged.size()<KEEP_FINAL)merged.push(c);else if(c.key<merged.top().key){merged.pop();merged.push(c);}}}
  while(!merged.empty()){auto c=merged.top();merged.pop();std::cout<<"22 "<<c.p<<" "<<c.q<<" "<<c.h<<" "<<double(c.score)/SCALE<<" "<<c.logN<<" "<<mode<<"\n";}
 }
 std::cerr<<"FINISHED_GRID "<<counterPairs<<" SECONDS "<<omp_get_wtime()-started<<"\n";
}
