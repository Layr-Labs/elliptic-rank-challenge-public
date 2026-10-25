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
 std::ifstream input("submission/research/z4-tables.txt");
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
 double three[]={0, 1.459237, 2.047788, 3.200999, 4.451785, 5.71214, 6.375587, 8.339871, 10.189149, 12.294657, 15.3894, 18.746574, 23.595951, 32.024466};
 for(int f=1;f<=2;f++)for(int r=1;r<=13;r++)baselines[f][r]=three[r];
 std::priority_queue<Candidate> best[7][14];
 long long examined=0;
 for(int f=1;f<=2;f++){
  int n=limit;
  int low=-n;
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
    
    examined++;
    double s=sums[p+n];
    if(s>-10||std::gcd(p,q)!=1)continue;
    long long aa,bb;
    if(f==1){aa=(8LL*p-q)*(32LL*p+7*q);bb=8LL*q*(p+q)*(15LL*p-8*q)*(31LL*p-7*q);}
    else{aa=-8LL*q*(80LL*p+9*q);bb=2LL*(p-2*q)*(2LL*p-q)*(18LL*p-q)*(2LL*p-81*q);}
    if(!aa||!bb||(__int128)aa*aa==16*(__int128)bb)continue;
    auto val=[](long long value,int prime){if(value==0)return 60;int n=0;while(value%prime==0){value/=prime;n++;}return n;};
    double reduction=0;for(int prime:{2,3,5,7,11,13,19})reduction+=std::min(val(aa,prime),val(bb,prime)/2)*std::log(prime);
    long double a1=aa,a2=bb,a3=a1*a2,b2=a1*a1+4*a2,b4=a1*a3,b6=a3*a3;
    long double c4=b2*b2-24*b4,c6=-b2*b2*b2+36*b2*b4-216*b6;
    double h=std::max(logl(fabsl(c4))/4,logl(fabsl(c6))/6)-reduction;
    int r=1;while(r<13&&h>=baselines[f][r])r++;
    if(h>=20.0)continue;
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
