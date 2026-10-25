read("lib.gp");
PR3=primes([3,30000]);
\\ L=[[p,q,ln,H],...] -> p q S3 sc3
s3(L,out)={ for(i=1,#L, my(F=iferr(x400afamily(L[i][1],L[i][2]),e,[[]]), M, s=0.); if(#F[1]==0,next); M=ellminimalmodel(F[1]);
   for(j=1,#PR3, my(p=PR3[j]); s-=ellap(M,p)*log(p)/p);
   write(out, strprintf("%d %d %.3f %.3f", L[i][1], L[i][2], s, s-0.5*L[i][3]+L[i][4]))); }
