# Own find, lane z2-exh (ec-rank s1). Method: exhaustive enumeration of all Z/2 curves
# y^2 = x^3 + a x^2 + b x (reduced (a,b), both u=1 and u=2 minimal models) with naive height H in [3.797737, 4.852134),
# 2.47e11 curves, Mestre-Nagao sieve S1 = sum_{3<=p<=60} -a_p log p/p >= 6.5 then S2 (p<1000) >= 26 (35382 survivors),
# then PARI 2.17.4 ellrank on the minimal model. This curve: (a,b) = (-19695, 86618064), S2 rank 7 by ellrank.
# Code and logs: challenges/ec-rank/work/search/z2-exh/ (sieve.c, srun.sh, s7.txt, r7.out.*).
E 1 -1 1 -2667476 43710086
P -1628 9153
P -1514 25458
P -1290 37207
P 10 4122
P 1788 30538
P 3312 164269
P 4963 327718
#CLAIM rank=7
