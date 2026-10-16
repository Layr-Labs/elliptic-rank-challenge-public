# Own find, lane z2-r9 (ec-rank s1). Method: residue-class restricted enumeration of Z/2 curves
# y^2 = x^3 + a x^2 + b x (reduced (a,b), u=1 and u=2 minimal models) with naive height H in [0.0, 6.430661):
# b mod 3*5*7*11*13*17*19*23 classes with S_E = sum_{3<=p<=23} -a_p log p/p in [6.163, 1e+30) enumerated by CRT,
# then Mestre-Nagao S1 (p<60) >= 10.0, S2 (p<1000) >= 32.0, then PARI 2.17.4 ellrank on the minimal model.
# This curve: (a,b) = (-52867, -33936708768), S1 12.566 S2 39.431, ellrank rank >= 9.
# Code and logs: challenges/ec-rank/work/search/z2-r9/ (crt.c, run.py, s/p1.*).
E 1 1 0 -2179271791 -9516362821775
P -12214 3914987
P -24273/4 14957927/8
P 71844 14273153
P 131844 44599553
P -261165200/13689 8042049153035/1601613
P -372104/9 85565621/27
P 108127051/2116 342438086647/97336
P 730158841/14400 5515923501679/1728000
P 1071853037/16129 24094949245337/2048383
#CLAIM rank=9
