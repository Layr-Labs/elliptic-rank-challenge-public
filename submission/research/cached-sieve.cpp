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
struct Table {int prime; std::vector<double> weight; std::vector<int> inv;};
struct Candidate {int f,p,q,r; double h,s; bool operator<(const Candidate& b)const{return s<b.s;}};
int main(int argc,char**argv){
 int limit=argc>1?std::stoi(argv[1]):3000;
 std::vector<Table> tables[7];
 std::ifstream input("submission/research/modular-tables.txt");
 std::string line;
 while(std::getline(input,line)){
  for(char& c:line) if(c=='['||c==']'||c==',')c=' ';
  std::istringstream row(line);int f,p,a;row>>f>>p;
  Table t; t.prime=p;
  while(row>>a)t.weight.push_back(a*std::log(p)/p);
  if(t.weight.size()!=p)return 2;
  t.inv.resize(p);t.inv[1]=1;for(int j=2;j<p;j++)t.inv[j]=p-(p/j)*t.inv[p%j]%p;
  tables[f].push_back(t);
 }
 double baselines[7][8]={{}};
 double five[]={0,1.534971,2.893369,4.109783,5.840106,8.051843,9.988902,13.428116};
 double six[]={0,1.838270,2.919120,4.360639,5.446254,7.724108,9.959559,12.699047};
 for(int r=1;r<=7;r++){baselines[5][r]=five[r];baselines[6][r]=six[r];}
 std::priority_queue<Candidate> best[7][8];
 long long examined=0;
 for(int f=5;f<=6;f++){
  int n=f==5?limit:limit*2/3;
  int low=f==5?1:-n;
  for(int q=1;q<=n;q++){
   std::vector<double> sums(2*n+1,0.);
   for(const Table&t:tables[f]){
    int prime=t.prime;if(prime>251)continue;int qr=q%prime;if(!qr)continue;
    int inv=t.inv[qr];
    int residue=((low%prime+prime)%prime)*inv%prime;
    for(int p=low;p<=n;p++){
     sums[p+n]+=t.weight[residue];
     residue+=inv;if(residue>=prime)residue-=prime;
    }
   }
   for(int p=low;p<=n;p++){
    if(p==0)continue;
    if(f==6&&p==-q)continue;
    examined++;
    double s=sums[p+n];
    if(s>-10||std::gcd(p,q)!=1)continue;
    long double a1=q-p,a2=f==5?-(long double)p*q:-(long double)p*(p+q);
    long double a3=f==5?-(long double)p*q*q:-(long double)p*q*(p+q);
    long double b2=a1*a1+4*a2,b4=a1*a3,b6=a3*a3;
    long double c4=b2*b2-24*b4,c6=-b2*b2*b2+36*b2*b4-216*b6;
    double h=std::max(logl(fabsl(c4))/4,logl(fabsl(c6))/6);
    int r=1;while(r<7&&h>=baselines[f][r])r++;
    if(r==7&&h>=(f==5?13.199808:12.699047))continue;
    auto& heap=best[f][r];
    Candidate c={f,p,q,r,h,s};
    if(heap.size()<50000)heap.push(c);
    else if(s<heap.top().s){heap.pop();heap.push(c);}
   }
  }
  std::cerr<<"FAMILY "<<f<<" EXAMINED "<<examined<<"\n";
 }
 std::cout<<std::setprecision(17);
 for(int f=5;f<=6;f++)for(int r=1;r<=7;r++){
  std::priority_queue<Candidate> refined;
  auto& heap=best[f][r];
  while(!heap.empty()){
   auto c=heap.top();heap.pop();double score=0;
   for(const Table&t:tables[f]){
    int prime=t.prime,qr=c.q%prime;if(!qr)continue;
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
