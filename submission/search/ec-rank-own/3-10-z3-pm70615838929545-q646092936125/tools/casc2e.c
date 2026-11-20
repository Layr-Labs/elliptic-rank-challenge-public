/* casc2e.c -- casc2 with early exit for -t T (Theta_1 rows one Selmer element at a time).
 * casc2.c -- second descent by 2-isogeny (Fisher 2017, level 1: the alternating pairing Theta_1 on each 2-isogeny
 * Selmer group) for z10, in C on libpari.  Same construction as cas2.gp (validated against PARI ellrank r2):
 *   E: y^2 = x(x^2 + A x + B), E': y^2 = x(x^2 - 2A x + A^2 - 4B);  for xi in the phi-Selmer group of (A,B):
 *   conic xi r^2 = xi^2 S^2 + A xi S T + B T^2 (qfsolve, primes supplied), (S : T : r) = (f : g : h)(l, m) (qfparam, reduced),
 *   C_xi = {y^2 = f g}, pushout f(l,m)/m^2;  Theta_1(xi, eta) = sum_v (f(l_v), eta)_v with f g (l_v) a local square.
 *   S2 = kernel; bound = dim S2(E) + dim S2(E') - 2.  Selmer bases from selmer_core.h (bphi_z10_ex, validated).
 * Build: gcc -O2 -march=native -I../../work/selmer/pari/include -o casc2 casc2.c ../../work/selmer/pari/libpari.a -lm
 * Usage: casc2 [-t T] [-b] < 'fam p q ...' (fam z10; z10b = 5-isogenous pair (E5,E5'), also -b on z10 lines; z8 = pair E0 / E0/<4P>; z12 = pair E / E/<6P>)
 *        ->  'fam p q s1 s2 s1p s2p bound ms'   (-t: keep lines with bound >= T, + ' C<bound>')
 */
#include <pari/pari.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#undef ulong
#include "selmer_core.h"
#include "tables.h"

static u32 *spf; static u64 spf_lim;
static seltab_t TT; static u32 qroff_h[QRMAX]; static u64 *qrbits_h;
static u64 rs_ = 0x9E3779B97F4A7C15ULL;
static inline u64 rnd64(void) { rs_ ^= rs_ << 13; rs_ ^= rs_ >> 7; rs_ ^= rs_ << 17; return rs_; }
static inline double nowms(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return 1e3 * t.tv_sec + 1e-6 * t.tv_nsec; }
static double T_conic = 0, T_loc = 0;
#define NUMQ(t) (typ(t) == t_INT ? (t) : gel(t, 1))
#define DENQ(t) (typ(t) == t_INT ? gen_1 : gel(t, 2))

static GEN mask_gen(u64 m, const u64 *B, int n) {
    GEN r = gen_1;
    for (int q = 1; q < n; q++) if ((m >> q) & 1) r = mulii(r, utoipos(B[q]));
    if (m & 1) r = negi(r);
    return r;
}
/* F2 kernel of rows R (n columns) -> masks */
static int kernel_masks(const u64 *R, int nr, int n, u64 *out) {
    u64 piv[64]; int pc[64]; int rk = 0; u64 pm = 0;
    for (int i = 0; i < nr; i++) { u64 y = R[i];
        while (y) { int c = __builtin_ctzll(y); if (!((pm >> c) & 1)) { piv[rk] = y; pc[rk] = c; pm |= 1ULL << c; rk++; break; }
            for (int k = 0; k < rk; k++) if (pc[k] == c) { y ^= piv[k]; break; } } }
    /* reduce to RREF */
    for (int k = 0; k < rk; k++) for (int k2 = 0; k2 < rk; k2++) if (k2 != k && ((piv[k2] >> pc[k]) & 1)) piv[k2] ^= piv[k];
    int nk = 0;
    for (int f = 0; f < n; f++) if (!((pm >> f) & 1)) {
        u64 x = 1ULL << f;
        for (int k = 0; k < rk; k++) if ((piv[k] >> f) & 1) x |= 1ULL << pc[k];
        out[nk++] = x;
    }
    return nk;
}
static u64 mm(u64 a, u64 b, u64 m) { return (u64)((unsigned __int128)a * b % m); }

/* class of a nonzero rational F in (Q_v^x)/(Q_v^x)^2: alpha = v(F) mod 2, u = (v odd) non-residue flag of the unit, (v = 2) unit mod 8 */
static void qclass(GEN F, ulong v, int *alpha, int *u) {
    GEN pv = utoipos(v); long e = gvaluation(F, pv);
    GEN un = gdiv(F, powis(pv, e));
    ulong r = umodiu(NUMQ(un), v == 2 ? 8 : v), d = umodiu(DENQ(un), v == 2 ? 8 : v);
    *alpha = (int)(e & 1);
    if (v == 2) *u = (int)((r * d) & 7);                          /* d odd: d^-1 = d mod 8 */
    else { u64 w = mm(r, (u64)Fl_inv(d, v), v); *u = legendre(w, v) == -1; }
}
static int qsquare(GEN F, ulong v) {                               /* nonzero rational F a square in Q_v */
    int a, u; qclass(F, v, &a, &u);
    if (a) return 0;
    return v == 2 ? (u == 1) : !u;
}

/* exact v-adic point search (odd v) on y^2 = f g (l): disk tree l = a + b t.  In a disk, H(t) = fg(a + b t) = v^e G(t) with
 * G primitive; a residue t0 with G(t0) != 0 mod v gives points iff e is even and G(t0) is a square mod v (Hensel); residues
 * with G(t0) = 0 mod v are refined.  Small v: every residue; large v: random residues plus the roots of G mod v.
 * Returns the rational l or NULL.  The chart at infinity is handled by the caller (reversed polynomials). */
static GEN hensel_pt(GEN fg, GEN f, GEN a, GEN b, GEN pv, ulong v, int depth) {
    if (depth > 40) return NULL;
    GEN H = gsubst(fg, 0, gadd(a, gmul(b, pol_x(0))));
    long e = LONG_MAX; for (long i = 2; i < lg(H); i++) if (!gequal0(gel(H, i))) { long w = gvaluation(gel(H, i), pv); if (w < e) e = w; }
    if (e == LONG_MAX) return NULL;
    GEN G = gdiv(H, powis(pv, e));
    GEN Gb = cgetg(lg(G), t_POL); Gb[1] = G[1];
    for (long i = 2; i < lg(G); i++) gel(Gb, i) = Rg_to_Fp(gel(G, i), pv);
    Gb = FpX_renormalize(Gb, lg(Gb));
    int full = v < 2000;
    long ntry = full ? (long)v : 256;
    if (!(e & 1)) for (long k = 0; k < ntry; k++) {
        ulong t0 = full ? (ulong)k : rnd64() % v;
        GEN gv = FpX_eval(Gb, utoi(t0), pv);
        if (!signe(gv) || kronecker(gv, pv) != 1) continue;
        for (int sh = 0; sh < 3; sh++) {   /* f(x) = 0 is avoided by moving inside the disk */
            GEN x = gadd(a, gmul(b, addii(utoi(t0), mului(sh, pv))));
            if (!gequal0(poleval(f, x))) return x;
        }
    }
    if (degpol(Gb) <= 0) return NULL;
    GEN R = FpX_roots(Gb, pv);
    for (long i = 1; i < lg(R); i++) {
        GEN x = hensel_pt(fg, f, gadd(a, gmul(b, gel(R, i))), gmul(b, pv), pv, v, depth + 1);
        if (x) return x;
    }
    return NULL;
}
/* local value of the pushout function at v: returns class via *alpha, *u (v = 0: *u = sign bit) */
static int locval(GEN f, GEN g, ulong v, int *alpha, int *u) {
    pari_sp av = avma;
    if (v == 0) {
        GEN fg = gmul(f, g); long prec = nbits2prec(256 + 4 * gexpo(fg));
        GEN rr = realroots(fg, NULL, prec);
        GEN cand = cgetg(1, t_VEC);
        long nrr = lg(rr) - 1;
        if (nrr == 0) cand = mkvec3(gen_0, gen_1, gen_m1);
        else {
            GEN M = gaddsg(1, gsupnorm(rr, prec));
            cand = vec_append(cand, gmulsg(-2, M)); cand = vec_append(cand, gmulsg(2, M));
            for (long i = 1; i < nrr; i++) cand = vec_append(cand, gmul2n(gadd(gel(rr, i), gel(rr, i + 1)), -1));
        }
        for (long i = 1; i < lg(cand); i++) {
            GEN x = bestappr(gel(cand, i), int2n(prec2nbits(prec) / 2));
            GEN F = poleval(f, x), FG = poleval(fg, x);
            if (!gequal0(F) && gsigne(FG) > 0) { *alpha = 0; *u = gsigne(F) < 0; set_avma(av); return 0; }
        }
        set_avma(av); return -1;
    }
    if (v > 2) {   /* fast path: l mod v with f g (l) a nonzero square mod v */
        GEN fc = RgX_to_RgC(f, 3), gc = RgX_to_RgC(g, 3);
        u64 f0 = umodiu(gel(fc, 1), v), f1 = umodiu(gel(fc, 2), v), f2 = umodiu(gel(fc, 3), v);
        u64 g0 = umodiu(gel(gc, 1), v), g1 = umodiu(gel(gc, 2), v), g2 = umodiu(gel(gc, 3), v);
        for (int t = 0; t < 128; t++) {
            u64 l = rnd64() % v;
            u64 fv = (f0 + mm(f1, l, v) + mm(f2, mm(l, l, v), v)) % v, gv = (g0 + mm(g1, l, v) + mm(g2, mm(l, l, v), v)) % v;
            if (!fv || !gv) continue;
            if (legendre(mm(fv, gv, v), v) != 1) continue;
            *alpha = 0; *u = legendre(fv, v) == -1; set_avma(av); return 0;
        }
    }
    if (v > 2) {   /* exact disk tree, affine chart then the chart at infinity (l = 1/t: f g (1/t) t^4, f(1/t) t^2) */
        GEN fg = gmul(f, g), pv = utoipos(v);
        GEN x = hensel_pt(fg, f, gen_0, gen_1, pv, v, 0);
        if (x) { qclass(poleval(f, x), v, alpha, u); set_avma(av); return 0; }
        GEN fc = RgX_to_RgC(f, 3), gc = RgX_to_RgC(g, 3);
        GEN fr = deg2pol_shallow(gel(fc, 1), gel(fc, 2), gel(fc, 3), 0), gr = deg2pol_shallow(gel(gc, 1), gel(gc, 2), gel(gc, 3), 0);
        GEN fgr = gmul(fr, gr);
        x = hensel_pt(fgr, fr, gen_0, pv, pv, v, 1);   /* t in v Z_v: l = 1/t outside Z_v */
        if (x) { qclass(poleval(fr, x), v, alpha, u); set_avma(av); return 0; }
    }
    /* rational sampling (v = 2, or the exact search failed) */
    GEN fg = gmul(f, g), pv = utoipos(v);
    long K = 2 * gvaluation(content(fg), pv) + 20;
    GEN vK = powuu(v, K);
    pari_sp av2 = avma;
    for (long tries = 0; tries < 200000; tries++, set_avma(av2)) {
        int m = (int)(rnd64() % 3); GEN x;
        if (m == 0) x = randomi(vK);
        else if (m == 1) x = gdiv(randomi(vK), powuu(v, rnd64() % 6));
        else x = mulii(powuu(v, rnd64() % K), addiu(randomi(powuu(v, 4)), 1));
        GEN F = poleval(f, x), FG = poleval(fg, x);
        if (gequal0(F) || gequal0(FG)) continue;
        if (!qsquare(FG, v)) continue;
        qclass(F, v, alpha, u); set_avma(av); return 0;
    }
    set_avma(av); return -1;
}

/* Theta_1 rank on the Selmer group given by kernel masks km[0..nb) of (A,B) over basis B_[0..n); -1 on error */
/* c0, w: disc A^2 - 4B = c0 w^2 with c0 built from the forms (all its primes are basis primes) */
/* maxrk >= 0: early exit as soon as the rank of the rows known so far exceeds maxrk (*early = 1; returns that lower bound) */
static long n_rows_done = 0, n_rows_full = 0, n_early = 0, n_diag = 0;
static int g_diag = 0;   /* set when a computed Theta row has a nonzero diagonal entry (review R1): no early exit */
static int theta1_rank(GEN A, GEN Bc, GEN c0, GEN w, const u64 *B_, int n, const u64 *km, int nb, int maxrk, int *early) {
    pari_sp av = avma;
    *early = 0;
    GEN Pl = cgetg(1, t_VEC); for (int q = 1; q < n; q++) Pl = vec_append(Pl, utoipos(B_[q]));
    u64 H[64];
    int ord[64]; double cost[64];
    for (int a = 0; a < nb; a++) { double c = 0; for (int q = 1; q < n; q++) if ((km[a] >> q) & 1) c += log2((double)B_[q]); cost[a] = c; ord[a] = a; }
    for (int i = 1; i < nb; i++) { int x = ord[i], j = i; while (j > 0 && cost[ord[j - 1]] > cost[x]) { ord[j] = ord[j - 1]; j--; } ord[j] = x; }
    u64 epiv[64]; u64 epm = 0; int erank = 0, altok = 1;
    n_rows_full += nb;
    for (int ia = 0; ia < nb; ia++) {
        int a = ord[ia];
        pari_sp av2 = avma;
        double t0 = nowms();
        GEN xi = mask_gen(km[a], B_, n);
        GEN G = mkmat3(mkcol3(shifti(sqri(xi), 1), mulii(A, xi), gen_0), mkcol3(mulii(A, xi), shifti(Bc, 1), gen_0), mkcol3(gen_0, gen_0, shifti(negi(xi), 1)));
        /* small diagonal conic X^2 - c0 Y^2 - xi Z^2 = 0  (X = 2 xi^2 S + A xi T, Y = w xi T, Z = 2 xi r) */
        GEN Gs = mkmat3(mkcol3(gen_1, gen_0, gen_0), mkcol3(gen_0, negi(c0), gen_0), mkcol3(gen_0, gen_0, negi(xi)));
        GEN ss = qfsolve(mkvec2(Gs, Pl));
        if (typ(ss) != t_COL) { set_avma(av); return -1; }
        GEN X = gel(ss, 1), Y = gel(ss, 2), Z = gel(ss, 3);
        GEN Tt = gdiv(Y, mulii(w, xi)), rr = gdiv(Z, shifti(xi, 1));
        GEN Ss = gdiv(gsub(X, gmul(mulii(A, xi), Tt)), shifti(sqri(xi), 1));
        GEN sol = Q_primpart(mkcol3(Ss, Tt, rr));
        if (!gequal0(qfeval(G, sol))) { set_avma(av); return -4; }
        GEN M = qfparam(G, sol, 1);
        GEN f = deg2pol_shallow(gcoeff(M, 1, 1), gcoeff(M, 1, 2), gcoeff(M, 1, 3), 0);
        GEN g = deg2pol_shallow(gcoeff(M, 2, 1), gcoeff(M, 2, 2), gcoeff(M, 2, 3), 0);
        T_conic += nowms() - t0; t0 = nowms();
        /* places: infinity, the basis primes, and primes of the contents of f and g */
        GEN extra = gtoset(shallowconcat(gel(Z_factor(absi(content(f))), 1), gel(Z_factor(absi(content(g))), 1)));
        u64 h = 0; int al, u;
        if (locval(f, g, 0, &al, &u)) { set_avma(av); return -2; }
        if (u) h ^= 1;
        long nex = lg(extra) - 1;
        for (long j = 1; j < n + nex; j++) {
            ulong v; int jb;
            if (j < n) { v = (ulong)B_[j]; jb = j; }
            else { v = itou(gel(extra, j - n + 1)); jb = -1; int dup = 0; for (int q = 1; q < n; q++) if (B_[q] == v) dup = 1; if (dup) continue; }
            if (locval(f, g, v, &al, &u)) { set_avma(av); return -3; }
            u64 hm = 0;
            if (v == 2) {
                int e3 = (u & 3) == 3, w3 = (u == 3 || u == 5);
                if (e3) hm |= 1; if (w3) hm |= 2;
                for (int q = 2; q < n; q++) { u64 b = B_[q]; int be3 = (b & 3) == 3, bw = ((b & 7) == 3 || (b & 7) == 5); if ((e3 & be3) ^ (al & bw)) hm |= 1ULL << q; }
            } else {
                if (al) {   /* (F, b)_v = (b|v)^alpha for b != v */
                    if ((v & 3) == 3) hm |= 1;
                    { int r8 = (int)(v & 7); if (r8 == 3 || r8 == 5) hm |= 2; }
                    for (int q = 2; q < n; q++) if (q != jb && legendre(B_[q] % v, v) == -1) hm |= 1ULL << q;
                }
                if (jb >= 0) { int self = (al && (v & 3) == 3) ^ u; hm = (hm & ~(1ULL << jb)) | ((u64)self << jb); }
            }
            h ^= hm;
        }
        T_loc += nowms() - t0;
        H[a] = h;
        set_avma(av2);
        n_rows_done++;
        if (maxrk >= 0) {
            u64 r = 0; for (int b = 0; b < nb; b++) if (__builtin_parityll(H[a] & km[b])) r |= 1ULL << b;
            if ((r >> a) & 1) { altok = 0; if (!g_diag) n_diag++; g_diag = 1; }   /* nonzero self-pairing: bad row -> no rounding, no early exit */
            while (r) { int c = __builtin_ctzll(r); if (!((epm >> c) & 1)) { epiv[c] = r; epm |= 1ULL << c; erank++; break; } r ^= epiv[c]; }
            int lb = altok ? ((erank + 1) & ~1) : erank;   /* Theta_1 is alternating: its rank is even */
            if (!g_diag && lb > maxrk) { *early = 1; set_avma(av); return lb; }
        }
    }
    /* Theta matrix rows: bit b = parity(H[a] & km[b]); rank */
    u64 Cr[64];
    for (int a = 0; a < nb; a++) { u64 r = 0; for (int b = 0; b < nb; b++) if (__builtin_parityll(H[a] & km[b])) r |= 1ULL << b; Cr[a] = r; }
    int rk = 0;
    for (int col = 0; col < nb; col++) {
        int piv = -1; for (int r = rk; r < nb; r++) if ((Cr[r] >> col) & 1) { piv = r; break; }
        if (piv < 0) continue;
        u64 t = Cr[piv]; Cr[piv] = Cr[rk]; Cr[rk] = t;
        for (int r = 0; r < nb; r++) if (r != rk && ((Cr[r] >> col) & 1)) Cr[r] ^= Cr[rk];
        rk++;
    }
    set_avma(av);
    return rk;
}

#define FAM_Z10B 100   /* z10 second pair (E5, E5') = 5-isogenous to (E, E') */
typedef struct { int s1, s2, t1, t2, bound, err, early; } c2res_t;
static c2res_t casc2_fam(int fam, i64 p, i64 q, int target) {
    c2res_t R = {0, 0, 0, 0, 0, 0, 0};
    g_diag = 0;
    pari_sp av = avma;
    i64 F[MAXF]; int nf = family_forms(fam == FAM_Z10B ? FAM_Z10 : fam, p, q, F);
    if (nf <= 0) { R.err = 1; return R; }
    fact_t ff[MAXF];
    for (int m = 0; m < nf; m++) if (factor64(F[m] < 0 ? (u64)(-F[m]) : (u64)F[m], &ff[m], spf, spf_lim)) { R.err = 2; return R; }
    u64 B_[MAXB], R1[3 * MAXB + 8], R2[3 * MAXB + 8]; int n, n1, n2, sd[2];
    int b = (fam == FAM_Z10) ? bphi_z10_ex(p, q, ff, &TT, sd, B_, &n, R1, &n1, R2, &n2)
          : (fam == FAM_Z10B) ? bphi_z10b_ex(p, q, ff, &TT, sd, B_, &n, R1, &n1, R2, &n2)
                             : bphi_one_ex(fam, p, q, ff, &TT, sd, B_, &n, R1, &n1, R2, &n2);
    if (b < 0) { R.err = 3; return R; }
    u64 k1[64], k2[64]; int nb1 = kernel_masks(R1, n1, n, k1), nb2 = kernel_masks(R2, n2, n, k2);
    if (nb1 != sd[0] || nb2 != sd[1]) { R.err = 4; return R; }
    /* A, B exactly and the decompositions A^2 - 4B = c0a wa^2, A2^2 - 4 B2 = 16 B = c0b wb^2 */
    GEN P = stoi(p), Q = stoi(q), PQ = subii(P, Q), U = subii(shifti(P, 1), Q);
    GEN A, Bc, c0a, wa, c0b, wb;
    if (fam == FAM_Z10) {
        GEN f1 = addii(subii(shifti(sqri(P), 1), shifti(mulii(P, Q), 1)), sqri(Q));
        GEN k4 = subii(addii(addii(subii(shifti(powiu(P, 4), 2), mulsi(12, mulii(powiu(P, 3), Q))), mulsi(6, mulii(sqri(P), sqri(Q)))), shifti(mulii(P, powiu(Q, 3)), 1)), powiu(Q, 4));
        A = negi(mulii(f1, k4));
        GEN P2 = addii(subii(sqri(P), mulsi(3, mulii(P, Q))), sqri(Q));
        Bc = mulii(mulii(utoipos(16), mulii(powiu(P, 5), powiu(PQ, 5))), P2);
        GEN F4 = subii(subii(shifti(sqri(P), 2), shifti(mulii(P, Q), 1)), sqri(Q));
        c0a = mulii(mulii(Q, U), F4); wa = sqri(mulii(Q, U));
        c0b = mulii(mulii(P, PQ), P2); wb = shifti(sqri(mulii(P, PQ)), 4);
    } else if (fam == FAM_Z10B) {   /* E5: A5 = -f1 K4b, B5 = 16 p (p-q) P2^5;  B5' = q u F4^5 */
        GEN f1 = addii(subii(shifti(sqri(P), 1), shifti(mulii(P, Q), 1)), sqri(Q));
        GEN k4 = subii(addii(subii(addii(shifti(powiu(P, 4), 2), mulsi(36, mulii(powiu(P, 3), Q))), mulsi(66, mulii(sqri(P), sqri(Q)))), mulsi(26, mulii(P, powiu(Q, 3)))), powiu(Q, 4));
        A = negi(mulii(f1, k4));
        GEN P2 = addii(subii(sqri(P), mulsi(3, mulii(P, Q))), sqri(Q));
        Bc = mulii(mulii(utoipos(16), mulii(P, PQ)), powiu(P2, 5));
        GEN F4 = subii(subii(shifti(sqri(P), 2), shifti(mulii(P, Q), 1)), sqri(Q));
        c0a = mulii(mulii(Q, U), F4); wa = sqri(F4);
        c0b = mulii(mulii(P, PQ), P2); wb = shifti(sqri(P2), 4);
    } else if (fam == FAM_Z8) {   /* E0: A = K8, B = 16 w^4;  B' = q^2 u^4 Q */
        GEN W = mulii(P, PQ);
        A = addii(subii(addii(subii(shifti(powiu(P, 4), 3), shifti(mulii(Q, powiu(P, 3)), 4)), shifti(mulii(sqri(Q), sqri(P)), 4)), shifti(mulii(powiu(Q, 3), P), 3)), powiu(Q, 4));
        Bc = shifti(powiu(W, 4), 4);
        GEN Qf = addii(subii(shifti(sqri(P), 3), shifti(mulii(P, Q), 3)), sqri(Q));
        c0a = Qf; wa = mulii(Q, sqri(U));
        c0b = gen_1; wb = shifti(sqri(W), 4);
    } else {                      /* z12 E: A = A12, B = 16 p^6 (p-q)^6 A3^2;  B' = q^2 u^6 A2^3 A6 */
        static const long c12[9] = {24, -96, 216, -312, 288, -168, 60, -12, 1};
        A = gen_0; for (int k = 0; k <= 8; k++) A = addii(A, mulsi(c12[k], mulii(powiu(Q, k), powiu(P, 8 - k))));
        GEN A3 = addii(subii(mulsi(3, sqri(P)), mulsi(3, mulii(P, Q))), sqri(Q));
        GEN A2f = addii(subii(shifti(sqri(P), 1), shifti(mulii(P, Q), 1)), sqri(Q));
        GEN A6 = addii(subii(mulsi(6, sqri(P)), mulsi(6, mulii(P, Q))), sqri(Q));
        Bc = shifti(sqri(mulii(powiu(mulii(P, PQ), 3), A3)), 4);
        c0a = mulii(A2f, A6); wa = mulii(mulii(Q, powiu(U, 3)), A2f);
        c0b = gen_1; wb = shifti(mulii(powiu(mulii(P, PQ), 3), A3), 4);
    }
    GEN A2 = shifti(negi(A), 1), B2 = subii(sqri(A), shifti(Bc, 2));
    if (!equalii(mulii(c0a, sqri(wa)), B2) || !equalii(mulii(c0b, sqri(wb)), shifti(Bc, 4))) { R.err = 6; set_avma(av); return R; }
    /* early exit for -t T: bound = nb1 + nb2 - 2 - r1 - r2 < T as soon as r1 + r2 > nb1 + nb2 - 2 - T */
    int mx = (target != -1000) ? nb1 + nb2 - 2 - target : -1, e1 = 0, e2 = 0;
    if (target != -1000 && mx < 0) { R.s1 = nb1; R.t1 = nb2; R.bound = nb1 + nb2 - 2; R.early = 1; n_early++; set_avma(av); return R; }
    int r1 = theta1_rank(A, Bc, c0a, wa, B_, n, k1, nb1, mx, &e1);
    if (r1 < 0) { R.err = 5; set_avma(av); return R; }
    if (e1) { R.s1 = nb1; R.t1 = nb2; R.s2 = nb1 - r1; R.t2 = nb2; R.bound = nb1 + nb2 - 2 - r1; R.early = 1; n_early++; set_avma(av); return R; }
    int r2 = theta1_rank(A2, B2, c0b, wb, B_, n, k2, nb2, mx >= 0 ? mx - r1 : -1, &e2);
    if (r2 < 0) { R.err = 5; set_avma(av); return R; }
    if (e2) { R.s1 = nb1; R.t1 = nb2; R.s2 = nb1 - r1; R.t2 = nb2 - r2; R.bound = R.s2 + R.t2 - 2; R.early = 1; n_early++; set_avma(av); return R; }
    R.s1 = nb1; R.s2 = nb1 - r1; R.t1 = nb2; R.t2 = nb2 - r2; R.bound = R.s2 + R.t2 - 2;
    set_avma(av);
    return R;
}

int main(int argc, char **argv) {
    int target = -1000;
    int sec = 0;   /* -b: run the second z10 pair (E5, E5') on z10 lines (same as family name z10b) */
    for (int i = 1; i < argc; i++) { if (!strcmp(argv[i], "-t") && i + 1 < argc) target = atoi(argv[++i]); else if (!strcmp(argv[i], "-b")) sec = 1; }
    pari_init(1L << 27, 0);
    setrand(utoipos(12345));
    spf_lim = 1u << 24; spf = (u32 *)calloc(spf_lim + 1, 4);
    for (u64 i = 2; i <= spf_lim; i++) if (!spf[i]) { spf[i] = (u32)i; if (i * i <= spf_lim) for (u64 j = i * i; j <= spf_lim; j += i) if (!spf[j]) spf[j] = (u32)i; }
    qrbits_h = (u64 *)calloc(QR_WORDS, 8); qr_build(qroff_h, qrbits_h, QR_WORDS);
    TT.tab2 = &tab2[0][0]; TT.tab3 = &tab3[0][0]; TT.qroff = qroff_h; TT.qrbits = qrbits_h; TT.tabw2 = 0;
    char line[4096], fs[32]; long long p, q; long nin = 0, nout = 0, nerr = 0; double ttot = 0;
    while (fgets(line, sizeof line, stdin)) {
        if (sscanf(line, "%31s %lld %lld", fs, &p, &q) != 3) continue;
        int fam = !strcmp(fs, "z10") ? FAM_Z10 : !strcmp(fs, "z10b") ? FAM_Z10B : !strcmp(fs, "z8") ? FAM_Z8 : !strcmp(fs, "z12") ? FAM_Z12 : -1;
        if (fam < 0) continue;
        if (sec && fam == FAM_Z10) fam = FAM_Z10B;
        nin++;
        if (nin % 1000 == 0) { fprintf(stderr, "progress %ld %s %lld %lld\n", nin, fs, p, q); fflush(stderr); }
        double t0 = nowms(); c2res_t R = casc2_fam(fam, p, q, target); double ms = nowms() - t0; ttot += ms;
        size_t L = strlen(line); while (L && (line[L - 1] == '\n' || line[L - 1] == '\r')) line[--L] = 0;
        if (R.err) { nerr++; fprintf(stderr, "error %d on %s %lld %lld\n", R.err, fs, p, q); }
        if (target == -1000) printf("%s %lld %lld %d %d %d %d %d %.3f%s\n", fs, p, q, R.s1, R.s2, R.t1, R.t2, R.err ? 99 : R.bound, ms, R.err ? " ERR" : "");
        else if (R.err || R.bound >= target) { printf("%s C%d%s%s\n", line, R.err ? 99 : R.bound, R.err ? " ERR" : "", g_diag ? " DIAG" : ""); fflush(stdout); nout++; }
        if (g_diag) fprintf(stderr, "diag: nonzero self-pairing (bad row, full computation done) on %s %lld %lld, bound %d\n", fs, p, q, R.bound);
    }
    fprintf(stderr, "casc2: %ld pairs, %.3f ms/pair (conics %.3f, local %.3f), %ld errors%s\n", nin, nin ? ttot / nin : 0, nin ? T_conic / nin : 0, nin ? T_loc / nin : 0, nerr, "");
    if (target != -1000) fprintf(stderr, "casc2: kept %ld\n", nout);
    fprintf(stderr, "casc2: early exits %ld, Theta rows computed %ld of %ld, curves with a nonzero self-pairing (no early exit) %ld\n", n_early, n_rows_done, n_rows_full, n_diag);
    pari_close();
    return 0;
}
