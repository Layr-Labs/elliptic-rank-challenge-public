\\ Independent PARI/GP checks of a board curve file; the file name is passed in the environment variable ECFILE.
default(parisizemax, 2*10^9);
parse(fn)={
  my(L=readstr(fn), pts=List(), ainv=0, v);
  for(i=1, #L,
    if(#L[i] < 2, next);
    v = strsplit(L[i], " ");
    if(v[1]=="E", ainv = [eval(v[2]),eval(v[3]),eval(v[4]),eval(v[5]),eval(v[6])]);
    if(v[1]=="P", listput(pts, [eval(v[2]), eval(v[3])])));
  [ainv, Vec(pts)]
};
ecfile = getenv("ECFILE");
PP = parse(ecfile);
E = ellinit(PP[1]); pts = PP[2];
print("PARI/GP version: ", version());
print("ellinit ok; discriminant has ", #Str(abs(E.disc)), " digits; ellminimalmodel(E)==E: ", ellminimalmodel(E)==E);
T = elltors(E);
print("elltors = ", T);
R = ellrank(E, 1);
print("ellrank(E,1): lower bound ", R[1], ", upper bound ", R[2], ", parity-related field ", R[3], ", #points returned ", #R[4]);
print("ellrootno(E) = ", ellrootno(E));
S = ellsaturation(E, pts, 100);
print("ellsaturation(points, 100): ", #S, " points");
H = matrix(#S,#S,i,j,ellheight(E,S[i],S[j]));
print("regulator (det of the height matrix of the saturated points) = ", matdet(H));
g=0; forprime(p=2,1000, if(E.disc%p!=0, g++));
print("good reduction at ", g, " of 168 primes below 1000");
print("number of points in the file: ", #pts, "; all on the curve: ", vecmin(vector(#pts,i,ellisoncurve(E,pts[i]))));
