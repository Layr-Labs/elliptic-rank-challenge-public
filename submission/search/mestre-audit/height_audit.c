#define main msr_original_main
#include "../mestre-own/msr.c"
#undef main
int main(void){
HMAX=100; SMIN=-1e9;
int a[4][10]={{0,1,7,8,15,16,24,25,34,35},{0,2,8,10,20,22,30,32,40,42},{0,3,12,14,26,27,39,42,55,57},{0,4,8,12,20,28,36,40,48,53}};
for(int j=0;j<4;j++){for(int i=0;i<10;i++)a[j][i]*=4;eval_set(a[j]);}
}
