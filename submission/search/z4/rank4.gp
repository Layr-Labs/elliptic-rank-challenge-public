\\ input: u D s S1 S2 H [S3]; output: u D s S2 S3 H rlo rhi ainvs(minimal)
runrank4(fin,fout,Hmax)={my(F=fileopen(fin,"r"),l); while(l=filereadstr(F), my(x=apply(eval,strsplit(l," ")));
  my(u=x[1],D=x[2],s=x[3]); my(E=ellinit([D*s,-D*u,-D^2*s*u,0,0])); if(!E.disc,next);
  my(M=iferr(alarm(30,ellminimalmodel(E)),e,0)); if(!M,next);
  my(h=log(max(abs(M.c4)^3,M.c6^2))/12.); if(h>=Hmax, next);
  my(R=iferr(alarm(60,ellrank(M,1)),e,[-2,-2,0,[]]));
  write(fout, u," ",D," ",s," ",x[5]," ",if(#x>6,x[7],0)," ",h," ",R[1]," ",R[2]," ",M[1..5])); fileclose(F)};
