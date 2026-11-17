# Focused Mestre audit

No board changes and no scoring curve found. All probes used GP 2.17.4 and small, explicitly chosen parameter lists. The existing `mestre-own` sources are unchanged.

## Useful defects and fixes

1. `msr.c:34-38` requires the quintic square-root polynomial to have integer coefficients, although rational coefficients are valid. An independent 1024-set, seed-260817 sample with roots in 0..72 passed only 142 sets. For the explicitly chosen roots `[0,11,17,33,47,51,59,62,65,68]`, coefficient denominator LCM is 256. Scaling every root by 4 makes that LCM 1 and preserves the exact minimal model and j-invariant. `dyadic.gp` verifies both polynomial identities and equality of minimal models. This is a coordinate change of the same independently found curve, not a separate board candidate. A correct rational-arithmetic implementation would also avoid the filter.

2. Coordinate scaling must be accounted for before the raw-height cutoff and fixed valuation ceiling. `msr.c:58` rejects raw H > HMAX+6 before local reduction; `lscale` caps 6e*log(p) at 42. This is not invariant under scaling roots. For `[0,1,7,8,15,16,24,25,34,35]`, exact minimal H is 12.85436693992907. The unchanged C height estimate gives 12.8544 after root scaling by 4 but 15.6270 after scaling by 32. Do not blindly increase root scales while retaining these filters. Remove the raw pre-cutoff and compute c4/c6 valuations with exact large integers, or divide out the known coordinate scaling first. In particular, `lscale` is not a substitute for `ellminimalmodel` at 2 and 3.

3. `mesPts` recovers only `hyperellratpoints(r,bnd)` points and skips both signs at the chosen base x. It does not explicitly seed known sections or the rational points at quartic infinity. Rational parameters with large coordinate numerators/denominators can have no points returned at a small bound while their algebraic sections still give rank >=12. `seeded.gp` directly maps both signs of all twelve root sections, the finite point above the other base sign, the infinity points when available, and Nagao's extra section. Verify all such added points exactly before independence filtering.

4. `denominator(polynomial)` in GP returned 1 for a polynomial with rational coefficients in this probe. Check `lcm(apply(denominator,Vec(q)))` instead.

## Construction and specialization checks

`audit.gp` proves the exact twelve-root condition. For centered

    F(x) = x^6+a*x^4+b*x^3+c*x^2+d*x+e,
    p(x,t) = F(x-t)*F(x+t) = Q(x,t)^2-R(x,t),

one has

    coefficient(R,x^5) = 12*t^2*(a*b-2*d).

Thus d=a*b/2 gives a quartic, subject to nonsingularity. The ten-root probes returned nine independent points, as expected; generic rank is around nine, so targeting rank 15-17 requires six to eight extra points.

The published Nagao family uses roots `[148,116,104,57,25,0] +/- t`. Its normalized quartic is

    (9*t^2+211950)*x^4 + (-2700*t^2-63901710)*x^3
    + (-18*t^4+396150*t^2+6706476489)*x^2
    + (2700*t^4-29575350*t^2-284435346600)*x
    + 9*t^6-159200*t^4+891699592*t^2+4156297690000.

Extra section:

    x=(t+703)/15,
    y=(-224*t^3-844*t^2+900484*t+2161725)/75.

The raw square-root-construction quartic is `(1248*t)^2` times this normalized quartic, so multiply this section's y by `1248*t` when using `mesqr`. A rational-infinity base change is

    t=(23550-u^2)/(2*u).

At independently chosen `u=77`, `t=-83/154`, the direct seed test found 13 independent points, determinant ~1.045850043604951e13, but H=18.648792, far too large for the target cells. At t=1 the simple bounded probe found 12 points at H=9.961191, also below required rank and above the rank-15 target height. Source: https://da.lib.kobe-u.ac.jp/da/kernel/E0003610/E0003610.pdf .

The complete Elkies rank-17 family coefficients and sections are in `elkies17.gp`, copied as *family formulas* from https://arxiv.org/html/2608.25406 . The exact polynomial identity for its first section passed. Fresh t=0,1,2,3 give H=15.9024,16.6479,17.1203,17.4777; the 17 displayed sections become dependent there. Fresh t=7/11 gives 17 independent points (smallest Gram eigenvalue ~0.09655) at H=21.222339, too large. The compact rank-18 base-change equation is

    U^2=4225*t^2+38636*t+289444,
    t=(289444-v^2)/(130*v-38636), U=65*t+v.

The conic identity passes exactly. The paper does not print an explicit eighteenth point for this conic, so that section still needs reconstruction; these formulas alone are not a ready 18-point certificate.

## CPU recommendation

Do not spend the next short allocation on another unmodified ten-root broad scan. If retaining it, first fix rational square-root coverage, height filtering, and explicit seeds; then constrain the exact minimal H before expensive point search, use deeper prime sums only on the resulting shortlist, and seed `ellrank` with the known points. No empirically defensible success-time estimate for a rank-15 curve below H=7.68185 follows from this tiny audit. The current fixed-height torsion family with precomputed finite-field traces has a better immediate CPU allocation than the high-height Elkies family. No full novelty-enabled verifier was run because none of these probes meets a target cell.
