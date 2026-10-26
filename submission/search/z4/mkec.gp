read("desc.gp");
\\ mkec(a, r, fn, hdr): a = minimal ainvs; find r independent points; write .ec file
mkec(a,r,fn,hdr)={my(E=ellinit(a)); my(R=ellrank(E,2)); my(pts=R[4]);
  if(#pts<r, my(D=desc2(E,20000)); pts=indep(E,concat(pts,D)));
  if(#pts<r, print("only ",#pts," points"); return(0));
  pts=pts[1..r];
  \\ try to make points smaller via LLL on height lattice (ellsaturation keeps it simple)
  pts=iferr(ellsaturation(E,pts,100),e,pts);
  my(F=fileopen(fn,"w"));
  filewrite(F,hdr);
  filewrite(F,Str("E ",a[1]," ",a[2]," ",a[3]," ",a[4]," ",a[5]));
  foreach(pts,P, filewrite(F,Str("P ",P[1]," ",P[2])));
  filewrite(F,Str("#CLAIM rank=",r)); fileclose(F); 1};
