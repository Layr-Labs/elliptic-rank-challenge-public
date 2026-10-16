# Own find, lane z2-r9 (ec-rank s1). Method: residue-class restricted enumeration of Z/2 curves
# y^2 = x^3 + a x^2 + b x (reduced (a,b), u=1 and u=2 minimal models) with naive height H in [0.0, 7.363717):
# b mod 3*5*7*11*13*17*19*23 classes with S_E = sum_{3<=p<=23} -a_p log p/p in [5.8, 1e+30) enumerated by CRT,
# then Mestre-Nagao S1 (p<60) >= 10.0, S2 (p<1000) >= 32.0, then PARI 2.17.4 ellrank on the minimal model.
# This curve: (a,b) = (-960398, 213690447032), S1 12.202 S2 45.519, ellrank rank >= 10.
# Code and logs: challenges/ec-rank/work/search/z2-r9/ (crt.c built as crt_hc5, r10c/run10.py, s/r10c.*).
E 0 1 0 -93764325769 2791717197247671
P -176125 117654536
P -31694 75707647
P -4105 56360896
P 1970 51058879
P 5417 47790740
P 968279 905442640
P -8513095/49 40348678628/343
P -479021/9 2358679456/27
P 189980093/4 2618503314325/8
P 8768043119/529 820880681905064/12167
#CLAIM rank=10
