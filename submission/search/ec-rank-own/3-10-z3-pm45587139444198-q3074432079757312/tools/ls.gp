\\ Z/4 family [1,-t,-t,0,0], t = u/(D s^2): a_l of the global minimal model when l | s (l >= 5, l not | u D): find the rule.
{ setrand(3); my(ok = 0, bad = 0, ex = List());
  for(k = 1, 600, my(l = prime(3 + random(40)), D = [1,2,3,5,6,7,10,11,13,14,15,-1,-2,-3,-6,-7][random(16)+1], s = l * (random(10^4) + 1), u = random(10^9) + 1);
    if(gcd(u, D * s) != 1 || u % l == 0 || D % l == 0, next);
    my(E = ellinit([D*s, -D*u, -D^2*s*u, 0, 0]), a = ellap(E, l), r1 = kronecker(-u*D, l), r2 = kronecker(u*D, l), r3 = kronecker(u, l), r4 = kronecker(-D, l), r5 = kronecker(D, l), r6 = kronecker(-u, l));
    listput(ex, [a, r1, r2, r3, r4, r5, r6]));
  my(names = ["-uD", "uD", "u", "-D", "D", "-u"]);
  for(j = 2, 7, print(names[j-1], ": matches ", #select(v -> v[1] == v[j], Vec(ex)), " of ", #ex));
  print("a_l values seen: ", Set(apply(v -> v[1], Vec(ex)))); }
quit
