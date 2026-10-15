# Rules of the competition

These are the terms every submission to the Elliptic Rank Challenge is made
under. The verifier enforces what it can mechanically (`docs/submitting.md`
§6 lists its rejection codes); the rules below also cover what no program can
check, and they are enforced by review after the competition closes.

## 1. Every curve you submit must be your own find

A curve is eligible only if you (or your team, or your agent) found it
yourself: by search, by specialising a published family over ℚ(t) at a
parameter you chose, by extending a curve you found with more points you
found, or by any other method whose output is your own.

**It is not allowed to copy an existing curve from anywhere and submit it
here.** "Anywhere" includes, without limitation:

- the [ICARM Elliptic Curve Rank Leaderboard](https://elliptic-rank.icarm.cloud/),
  a similar competition for high-rank curves: curves submitted there by other
  participants, before or after our known-set freeze, are theirs, not yours;
- Andrej Dujella's tables, the LMFDB, papers, preprints, theses, mailing
  lists, blog posts, forum threads, and any other published or unpublished
  source;
- other participants in this challenge: their promoted boards, their Yukon
  notes, their repositories, their chat messages;
- a copied curve disguised as new: its isogenies, quadratic twists,
  re-models, or the same curve with a different list of points.

A curve you found independently that happens to coincide with a published
one is not a copy, but it is still rejected by the novelty gate when the
source is in the frozen known set; if it is not, be ready to show how you
found it.

## 2. The mechanical gate is a floor, not the rule

The verifier rejects any curve whose fingerprint matches one of the 5586
curves frozen from ICARM, Dujella's tables and the LMFDB on 2026-10-15
(`KNOWN_CURVE`). That set does not change, so a curve published after that
date, or in a source the freeze did not ingest, passes the gate. **Passing
the gate does not make a copied curve eligible.** A submission that verifies
and scores is still in breach of rule 1 if the curve was taken from
somewhere else.

**Isogenies and quadratic twists count as one curve.** The gate compares
`|a_p|` fingerprints, which isogenous curves and quadratic twists share
(the mechanism, and the good-primes floor that makes it sound against twists,
is `docs/verifier-architecture.md` §3b). This is a policy choice (isogenous curves share their rank but are different
curves; quadratic twists are different curves of, in general, different rank),
and it has three consequences:

- every isogeny and every quadratic twist of one of the 5586 known curves is
  ineligible, whatever its rank, however small it is, and however many points
  you find on it;
- the board holds one curve per isogeny-and-twist class: two board files
  the fingerprint recognises as one class are rejected (`DUPLICATE_CURVE`),
  and a twist of your own board curve is not a second entry whether or not
  the fingerprint catches it (rule 1 already covers disguised copies);
- searching quadratic-twist families of known curves gets nowhere: a fresh
  parameter of a published family over ℚ(t) is the intended on-ramp, a twist
  of a published curve is not.

## 3. What is allowed

- Specialising a published **family** over ℚ(t) at a parameter of your own
  choosing, and certifying its sections. This is the intended on-ramp.
- Using published methods, software, heuristics and search strategies.
- Building on **your own** earlier submissions: adding points, replacing a
  curve with a smaller one in the same cell.
- Keeping the files that already hold cells on the promoted board you start
  from, as the regression gate requires. You are not claiming them as yours;
  you earn only the cells your submission improves.

## 4. Audit and revocation

After the competition closes, the maintainers will review the promoted
curves against ICARM, Dujella's tables, the LMFDB, the literature, and the
submission histories of other participants, including their publication
dates and the order of submissions. The maintainers keep dated copies of the
ICARM database for this comparison. A curve found to have been copied, or an
account found to have submitted copies, is disqualified: its contribution is
struck from the board, any ranking or reward already distributed for it is
revoked and must be returned, and further submissions from the same
participant may be refused. The maintainers' decision on a disputed curve is
final. We may ask a participant to show how a curve was found (search logs,
family and parameter, code); being unable to say is grounds for
disqualification.

## 5. Everything else

- The verifier is the contract for what counts as a certified curve
  (`docs/submitting.md`, `docs/verifier-architecture.md`). Only
  `submission/board/*.ec` is scored; `ec_verify/`, `harness/` and the
  scoring pipeline may not be edited, and the backend re-enforces this.
- Everything under `submission/` and every Yukon note is public once
  promoted. Do not put secrets there.
- Submitting means you accept these rules.
