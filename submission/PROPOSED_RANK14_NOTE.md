# A new promoted-board cell with fourteen certified points

## Proposed contribution

This submission adds one independently generated elliptic curve with rational
torsion group Z/2Z x Z/2Z and fourteen certified independent points of infinite
order. Its minimal-model naive height is 27.043057. It would hold cell 2x2/14,
which has no holder on the current promoted board. The frozen baseline for
that cell is 28.082546.

The complete local board verifies at 136.43579073440776 against the current
promoted parent total 132.73424244886414. The verified improvement is
3.7015482855436233 points. Only cell 2x2/14 improves, from contribution zero
to 3.701548285543624. This is an improvement under the challenge's frozen
height table, rather than a claim of the largest known rational rank.

The new scored file is submission/board/2x2-14-p96576-q9845.ec.
All twenty-six inherited curves remain byte-for-byte unchanged. The inherited
rank-13 holder is preserved in its own cell. No inherited curve, point list,
or contribution is claimed as a new discovery in this submission.

## Fresh construction and attribution

The new curve is a fresh specialization of the published X400a family at
the independently searched parameter T=96576/9845. For coprime integers p
and q with q positive, the family equations used here are

    A=q(121q-99p)(26p-23q)(75p-44q)
    B=p(21q-74p)(25p-36q)(p+77q)
    y^2=x(x+A)(x+B)

Here p=96576 and q=9845. The raw a-invariants are [0,A+B,0,A*B,0].
Exact minimalization gives the a-invariants

    [1, 0, 1,
     -1983174073675387598123283689382830167181410751,
     33210886676652029621965122111543787426920462301469572647735614824898]

The family equations were supplied in the header of the inherited file
submission/board/2x2-12-own.ec from Yukon submission
43af8780-369f-484f-b8df-47147017654b, where they are attributed to wgxli.
This source is credited as the origin of the family method. Its particular
curve and point list were not used as this new candidate.

Our earlier submission 7e9f6673-cf5f-4cb9-82b4-a22add86d720 published four
independently derived section formulas and a seeded PARI/GP point-search
procedure for this family. This result reuses that method at a fresh
parameter. The published section source used for reproduction is pinned at

https://github.com/Layr-Labs/elliptic-rank-challenge-public/blob/7482bf99c4a605a688f26c12d413f7085233dec2/submission/research/x400a-family.gp

The current inherited parent is our promoted submission
e4015be6-ebb5-475f-b505-f1b867f7cc35, commit
c1890e67b75e1de012b04c7add6c860ec2ba495e. It provides the board preserved
for the regression gate, including the smaller rank-13 curve. These three
material predecessors are explicitly declared through Yukon dependencies.

No individual published curve was copied as the new result. The rational
parameter, model and point witnesses for this specialization were computed
in the solver's own search. Published family formulas and already public
section methods are used as construction tools, as allowed by the challenge
rules. The claimed novelty still passes the official frozen-ledger check.

## Reproducing the model

Calculations used Ubuntu under WSL, the repository Python environment and
its pinned PARI/GP 2.17.4. After the normal challenge setup, these statements
reproduce the raw curve, its minimal model, torsion group and naive height:

    p=96576;q=9845;
    A=q*(121*q-99*p)*(26*p-23*q)*(75*p-44*q);
    B=p*(21*q-74*p)*(25*p-36*q)*(p+77*q);
    Eraw=ellinit([0,A+B,0,A*B,0]);
    E=ellminimalmodel(Eraw);
    print([E.a1,E.a2,E.a3,E.a4,E.a6]);
    print(elltors(E)[2]);
    print(log(max(abs(E.c4)^3,E.c6^2))/12);

The verifier recomputes these quantities from the submitted integral model;
a textual rank or height claim in a comment is not accepted as evidence.

For point generation, read the already published x400a-family.gp file at
the pinned commit above. Its x400afamily(p,q) function returns the raw
model and four rational section points. Run the following in the pinned GP:

    setrand(20261025);
    F=x400afamily(96576,9845);
    E=ellminimalmodel(F[1],&v);
    P=ellchangepoint(F[2],v);
    R=ellrank(E,0,P);
    if(R[2]>=14 && #R[4]<14,R=ellrank(E,4,R[4]));
    print([R[1],R[2]]);
    print(#R[4]);
    for(i=1,#R[4],print(R[4][i]));

A fresh process returned bounds [14,14] and regenerated the same fourteen
rational witnesses. The model and all coordinate strings were compared
exactly with the saved candidate. The official verifier then independently
accepted those witnesses with novelty checking enabled.

## Reproducing the certificate

The challenge certificate consists of the minimal integral a-invariants
and fourteen exact rational point pairs in the scored file. Reproducing
the search order is unnecessary to validate the result. The standard
verification commands are

    python -m ec_verify submission/board/2x2-14-p96576-q9845.ec
    python -m ec_verify submission/board --parent /path/to/parent/board

For the second command, the parent is the unmodified set of .ec files
from commit c1890e67b75e1de012b04c7add6c860ec2ba495e.
The local full-board command exited zero, reported parent_gate=checked,
and returned a nonempty improvements list containing only cell 2x2/14.

All fourteen submitted points remain after torsion stripping. The smallest
height-pairing eigenvalue is 2.130404159532254.
Agreement between the two required height precisions is
7.137051761873386e-56. The regulator lower bound is
2.66924848024946e+18. These values are recomputed by the
verifier and are not adjustable search parameters.

The checks cover exact point membership, integral minimality under the
challenge contract, the rational torsion group, independent infinite-order
points, frozen-ledger novelty, board duplicates and preservation of the
parent's held cells. The recorded rank means a certified lower bound of
fourteen from fourteen independent witnesses. The score relies on this
certificate, regardless of the separate GP rank-bound calculation.

## Score, scope and submitted contents

For a known baseline, the contribution is
100*(baseline-H)/baseline. Substituting baseline 28.082546 and height
27.043057 gives 3.701548285543624. The parent has no holder in this cell,
so this contribution is the new gain. The larger complete board total
includes inherited contributions and must not be confused with the gain
attributable to this individual new curve.

The submitted archive contains the scored board and its existing README.
The required public note covers the winning construction, reproduction,
certificate and provenance. Experimental sources, candidate lists,
unsubmitted curves and broader research logs are excluded from the archive.
The exact model attribution is the user-confirmed GPT-6.1 Sol, with Codex
as the coding harness. No separately confirmed effort setting or direct
human coauthor is asserted.

This finite search establishes a new board improvement for this curve.
It does not establish global height optimality in the family, exhaustive
coverage of rational parameters, or a new overall rank record. Inherited
holders are retained to satisfy the regression gate, and the verifier
remains the authority for both acceptance and scoring.
