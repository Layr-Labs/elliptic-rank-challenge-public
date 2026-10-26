\\ input: a b S1 S2 H S3; output: a b S2 S3 H rlo rhi ainvs(minimal) torsion
runrank3(fin,fout,Hmax)={my(F=fileopen(fin,"r"),l); while(l=filereadstr(F), my(x=apply(eval,strsplit(l," ")));
  my(E=ellinit([x[1],0,x[2],0,0])); if(!E.disc,next);
  my(M=iferr(alarm(30,ellminimalmodel(E)),e,0)); if(!M,next);
  my(h=log(max(abs(M.c4)^3,M.c6^2))/12.); if(h>=Hmax, next);
  my(R=iferr(alarm(60,ellrank(M,1)),e,[-2,-2,0,[]]));
  write(fout, x[1]," ",x[2]," ",x[4]," ",x[6]," ",h," ",R[1]," ",R[2]," ",M[1..5])); fileclose(F)};
