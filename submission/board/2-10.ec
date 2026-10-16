# Own find, lane z2-r9 (ec-rank s1). Method: residue-class restricted enumeration of Z/2 curves
# y^2 = x^3 + a x^2 + b x (reduced (a,b), u=1 and u=2 minimal models) with naive height H in [0.0, 7.363717):
# b mod 3*5*7*11*13*17*19*23 classes with S_E = sum_{3<=p<=23} -a_p log p/p in [5.8, 1e+30) enumerated by CRT,
# then Mestre-Nagao S1 (p<60) >= 10.0, S2 (p<1000) >= 32.0, then PARI 2.17.4 ellrank on the minimal model.
# This curve: (a,b) = (-459010, 157429335058), S1 12.634 S2 44.046, ellrank rank >= 10.
# Code and logs: challenges/ec-rank/work/search/z2-r9/ (crt.c built as crt_hc5, r10c/run10.py, s/r10c.*).
E 0 -1 0 87199275025 16923561771699111
P -98490 85906431
P 106845 165710844
P -154274754/1849 7570744014729/79507
P -858675/49 42551499108/343
P 692774/25 17393073243/125
P 16128366/49 97791248421/343
P 66448055/49 556141864242/343
P 131259905/4096 36841116547839/262144
P 178687663/361 2920980842302/6859
P 7540835757/5329 670848604605324/389017
#CLAIM rank=10
