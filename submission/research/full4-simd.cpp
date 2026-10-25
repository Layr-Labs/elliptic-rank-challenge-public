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
#include <omp.h>
constexpr int SCALE=128, BLOCK=32768, KEEP_LOCAL=25000, KEEP_FINAL=8000;
struct Table {int p;std::vector<int16_t>w;std::vector<int>inv;};
struct Candidate{int p,q,r,score;double h;bool operator<(const Candidate&b)const{return score<b.score;}};
using Heap=std::priority_queue<Candidate>;
int main(int argc,char**argv){
 int N=argc>1?std::stoi(argv[1]):150000,threads=argc>2?std::stoi(argv[2]):8;
 std::ifstream data("/tmp/ec-rank-full4-simd-tables.txt");std::string line;std::vector<Table>tables;
 int lower=0,upper=0;
 while(std::getline(data,line)){
  for(char&c:line)if(c=='['||c==']'||c==',')c=' ';
  std::istringstream row(line);Table t;row>>t.p;int a;
  while(row>>a)t.w.push_back(std::lround(SCALE*(a==9999?-0.4:a*std::log(t.p)/t.p)));
  if(t.w.size()!=t.p)return 2;
  lower+=std::min<int>(-51,*std::min_element(t.w.begin(),t.w.end()));
  upper+=std::max<int>(-51,*std::max_element(t.w.begin(),t.w.end()));
  t.inv.resize(t.p);t.inv[1]=1;for(int i=2;i<t.p;i++)t.inv[i]=t.p-(t.p/i)*t.inv[t.p%i]%t.p;
  tables.push_back(std::move(t));
 }
 if(lower<-32768||upper>32767){std::cerr<<"COUNTER_RANGE_UNSAFE "<<lower<<" "<<upper<<"\n";return 3;}
 std::cerr<<"PRIMES "<<tables.size()<<" COUNTER_BOUNDS "<<lower<<" "<<upper<<"\n";
 std::vector<std::vector<Heap>>heaps(threads,std::vector<Heap>(8));
 double bases[]={0,2.266502,3.212224,4.388376,6.219139,7.727479,10.740452,12.948901};
 omp_set_num_threads(threads);double started=omp_get_wtime();
 #pragma omp parallel
 {
  int tid=omp_get_thread_num();std::vector<int16_t>counts(BLOCK),period(2001);
  #pragma omp for schedule(dynamic,16)
  for(int q=1;q<N;q++){
   for(int lo=q+1;lo<=N;lo+=BLOCK){
    int len=std::min(BLOCK,N-lo+1);std::fill(counts.begin(),counts.begin()+len,0);
    for(const Table&t:tables){
     int prime=t.p,qr=q%prime;
     if(!qr){
      #pragma omp simd
      for(int i=0;i<len;i++)counts[i]-=51;
      continue;
     }
     int inv=t.inv[qr],residue=(int)((1LL*lo*inv)%prime);
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
     int score=counts[j];if(score>-14*SCALE)continue;
     int p=lo+j;if(std::gcd(p,q)!=1)continue;
     long double p2=(long double)p*p,q2=(long double)q*q;
     long double c4=16*(p2*p2+14*p2*q2+q2*q2);
     long double c6=64*(p2+q2)*(p2*p2-34*p2*q2+q2*q2);
     double h=std::max(logl(c4)/4,logl(fabsl(c6))/6);
     if((p&1)&&(q&1)){h-=std::log(2.);if((p-q)%8==0||(p+q)%8==0)h-=std::log(2.);}
     double hr=std::round(h*1e6)/1e6;
     int r=1;while(r<7&&hr>=bases[r])r++;
     if(hr>=bases[7])continue;
     auto&heap=heaps[tid][r];Candidate c={p,q,r,score,h};
     if(heap.size()<KEEP_LOCAL)heap.push(c);
     else if(score<heap.top().score){heap.pop();heap.push(c);}
    }
   }
   if(q%20000==0){
    #pragma omp critical
    std::cerr<<"DENOM "<<q<<" SECONDS "<<omp_get_wtime()-started<<"\n";
   }
  }
 }
 std::cout<<std::setprecision(17);
 for(int r=1;r<=7;r++){
  Heap merged;
  for(auto&byrank:heaps){auto&h=byrank[r];while(!h.empty()){auto c=h.top();h.pop();if(merged.size()<KEEP_FINAL)merged.push(c);else if(c.score<merged.top().score){merged.pop();merged.push(c);}}}
  while(!merged.empty()){auto c=merged.top();merged.pop();std::cout<<"42 "<<c.p<<" "<<c.q<<" "<<c.h<<" "<<double(c.score)/SCALE<<"\n";}
 }
 std::cerr<<"FINISHED_GRID "<<1LL*N*(N-1)/2<<" SECONDS "<<omp_get_wtime()-started<<"\n";
}
