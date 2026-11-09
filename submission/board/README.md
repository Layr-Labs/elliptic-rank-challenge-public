# Board

Put one certified curve per `*.ec` file in this directory (format: `docs/submitting.md` §2).
Only regular files named `[A-Za-z0-9._-]+.ec` directly in this directory are scored; anything
else here is ignored. Each file is verified from scratch; a curve counts in the cell of its
torsion group and certified rank, and the board total rises when a file undercuts that cell's
baseline. Never remove a file that holds a cell: the verifier rejects the regression.

Every curve you add must be your own find. Copying one from the ICARM leaderboard, a table, a
paper or another participant is against the rules even when the frozen known-curve gate does
not catch it; promoted curves are audited after the competition closes (`RULES.md`).
