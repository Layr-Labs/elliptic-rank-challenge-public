\\ Owner tool, 18 Nov 2026: triage of late candidates (class-wide bound >= target).  Driver sets INFILE (lines
\\ "fam p q ..."), OUTFILE, TARGET.  Per curve: H, BSD budget b = ellbsd(E) = Omega*prod(c_p)/#T^2, long Mestre-Nagao
\\ sums S(1e4), S(1e5), S(1e6), S(3e6) (all primes; an isogeny invariant), root number, and the regulator predicted
\\ by BSD if the rank equals TARGET: Reg*Sha = Lstar/b, b taken as the MAXIMUM over the isogeny class (search that
\\ member), Lstar = 4e4, 1e5, 1e6, 1e7 at rank <= 4, 5, 6, 7 (measured on saturated bases: review/budget/REPORT.md),
\\ Reg = 0.5 h^r; the tallest generator is typically 1.6 times that geometric mean.
\\ Reading: rank r curves show S(3e6) near (r - 1/2) * 7 to 9; a predicted mean height above 45 (B=1e5), 55 (1e6)
\\ or 64 (1e7) is out of reach of the 2-covering search.
default(parisizemax, 6000000000);
read("search/families/fam.gp");
{
my(L = readstr(INFILE));
foreach(L, s, if(#s == 0 || Vecsmall(s)[1] == 35, next);
  my(v = select(x -> #x, strsplit(s, " ")), f = v[1], p = eval(v[2]), q = eval(v[3]), E = minmodel(f, p, q), a = ellan(E, 3*10^6), S = 0., r = vector(4), k = 1, lim = [10^4, 10^5, 10^6, 3*10^6]);
  forprime(l = 2, 3*10^6, while(l > lim[k], r[k] = S; k++); S -= a[l]*log(l)/l); r[4] = S;
  my(b = ellbsd(E), bm = b, im = 1, L = ellisomat(E)[1]);
  for(k = 2, #L, my(bk = ellbsd(ellminimalmodel(ellinit(L[k][1])))); if(bk > bm, bm = bk; im = k));
  my(Ls = if(TARGET <= 4, 4e4, TARGET == 5, 1e5, TARGET == 6, 1e6, 1e7), R = Ls/bm, hm = (R/0.5)^(1./TARGET));
  write(OUTFILE, f, " ", p, " ", q, " H ", Strprintf("%.3f", heightE(E)), " w ", ellrootno(E), " log10b ", Strprintf("%.2f", log(b)/log(10)), " classmax ", Strprintf("%.2f", log(bm)/log(10)), " member ", im, "/", #L, " S1e4 ", Strprintf("%.1f", r[1]), " S1e5 ", Strprintf("%.1f", r[2]), " S1e6 ", Strprintf("%.1f", r[3]), " S3e6 ", Strprintf("%.1f", r[4]), " | if rank ", TARGET, ": Reg*Sha ~ ", Strprintf("%.1e", R), ", mean generator height ~ ", Strprintf("%.0f", hm)));
}
quit
