#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <queue>
#include <sstream>
#include <string>
#include <vector>
struct Table {int prime; std::vector<double> weight; std::vector<int> inv; double infinity;};
struct Candidate {int f,p,q,r; double h,s; bool operator<(const Candidate& b)const{return s<b.s;}};
int main(int argc,char**argv){
 int limit=argc>1?std::stoi(argv[1]):3000;
 std::vector<Table> tables[7];
 std::ifstream input("submission/research/strong-tables.txt");
 std::string line;
 while(std::getline(input,line)){
  for(char& c:line) if(c=='['||c==']'||c==',')c=' ';
  std::istringstream row(line);int f,p,a;row>>f>>p;
  Table t; t.prime=p;
  while(row>>a)t.weight.push_back(a==9999?-0.4:a*std::log(p)/p);
  if(t.weight.size()!=p+1)return 2;t.infinity=t.weight.back();t.weight.pop_back();
  t.inv.resize(p);t.inv[1]=1;for(int j=2;j<p;j++)t.inv[j]=p-(p/j)*t.inv[p%j]%p;
  tables[f].push_back(t);
 }
 double baselines[7][14]={{}};
 double three[]={0,1.213008,1.748024,2.404652,3.488487,4.177828,4.998527,6.663341,8.064684,9.658763,12.669669,14.744110,17.864161,20.199607};
 for(int f=1;f<=2;f++)for(int r=1;r<=13;r++)baselines[f][r]=three[r];
 long long roots[7][3],aa[7],bb[7],kk[7];
 std::ifstream config("submission/research/strong-coefficients.txt");int idx;
 while(config>>idx)config>>roots[idx][0]>>roots[idx][1]>>roots[idx][2]>>aa[idx]>>bb[idx]>>kk[idx];
 std::priority_queue<Candidate> best[7][14];
 long long examined=0;
 for(int f=1;f<=2;f++){
  int n=limit;
  int low=1;
  for(int q=1;q<=n;q++){
   std::vector<double> sums(2*n+1,0.);
   for(const Table&t:tables[f]){
    int prime=t.prime;if(prime>1000)continue;int qr=q%prime;if(!qr){for(int pp=low;pp<=n;pp++)sums[pp+n]+=t.infinity;continue;}
    int inv=t.inv[qr];
    int residue=((low%prime+prime)%prime)*inv%prime;
    for(int p=low;p<=n;p++){
     sums[p+n]+=t.weight[residue];
     residue+=inv;if(residue>=prime)residue-=prime;
    }
   }
   for(int p=low;p<=n;p++){
    if(p==0)continue;
    if(p==roots[f][0]*q||p==roots[f][1]*q||p==roots[f][2]*q)continue;
    examined++;
    double s=sums[p+n];
    if(s>-10||std::gcd(p,q)!=1)continue;
    long long p2=1LL*p*p,q2=1LL*q*q,a1i=aa[f]*p2+bb[f]*q2;
    long long c0=p2-roots[f][0]*roots[f][0]*q2,c1=p2-roots[f][1]*roots[f][1]*q2,c2=p2-roots[f][2]*roots[f][2]*q2;
    auto v2=[](long long value){if(value==0)return 60;return __builtin_ctzll((unsigned long long)std::llabs(value));};
    int shift=std::min(v2(a1i),(v2(kk[f])+v2(c0)+v2(c1)+v2(c2))/3);
    long double a1=a1i,a3=(long double)kk[f]*c0*c1*c2;
    long double b2=a1*a1,b4=a1*a3,b6=a3*a3;
    long double c4=b2*b2-24*b4,c6=-b2*b2*b2+36*b2*b4-216*b6;
    double h=std::max(logl(fabsl(c4))/4,logl(fabsl(c6))/6)-shift*std::log(2.);
    int r=1;while(r<13&&h>=baselines[f][r])r++;
    if(h>=20.199607)continue;
    auto& heap=best[f][r];
    Candidate c={f,p,q,r,h,s};
    if(heap.size()<50000)heap.push(c);
    else if(s<heap.top().s){heap.pop();heap.push(c);}
   }
  }
  std::cerr<<"FAMILY "<<f<<" EXAMINED "<<examined<<"\n";
 }
 std::cout<<std::setprecision(17);
 for(int f=1;f<=2;f++)for(int r=1;r<=13;r++){
  std::priority_queue<Candidate> refined;
  auto& heap=best[f][r];
  while(!heap.empty()){
   auto c=heap.top();heap.pop();double score=0;
   for(const Table&t:tables[f]){
    int prime=t.prime,qr=c.q%prime;if(!qr){score+=t.infinity;continue;}
    int residue=((c.p%prime+prime)%prime)*t.inv[qr]%prime;
    score+=t.weight[residue];
   }
   c.s=score;
   if(refined.size()<2500)refined.push(c);
   else if(score<refined.top().s){refined.pop();refined.push(c);}
  }
  while(!refined.empty()){
   auto c=refined.top();refined.pop();
   std::cout<<c.f<<" "<<c.p<<" "<<c.q<<" "<<c.h<<" "<<c.s<<"\n";
  }
 }
 std::cerr<<"TOTAL_GRID_PAIRS "<<examined<<"\n";
}
