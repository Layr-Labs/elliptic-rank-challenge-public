/* casc.c -- Cassels-Tate pairing on Sel2 for the full-2-torsion curves of z8 (E1) and z2z6, in C on libpari.
 * Same formula as cassels.gp (validated against PARI ellrank): for Lambda = (l1,l2,l3) in Sel2 (x - e_i = l_i u_i^2),
 * H_i: l_j u_j^2 - l_k u_k^2 = (e_k - e_j) u_0^2, q_i in H_i(Q) (qfsolve, primes supplied), L_i tangent at q_i;
 * <Lambda, M> = sum_v sum_i [ (L_i(x_v), m_i)_v = -1 ] mod 2, x_v in D_Lambda(Q_v), v | 2 r s (s-r) or v = inf.
 * Selmer basis from selmer_core.h (validated Legendre-symbol 2-descent).
 *
 * Build: gcc -O2 -march=native -I../../work/selmer/pari/include -o casc casc.c ../../work/selmer/pari/libpari.a -lm
 * Usage: casc [-t T] < lines 'fam p q ...' (fam z8 = E1, z2z6, z12, z2z8 = its 3 full-2-torsion curves)
 *        -> 'fam p q dim s r2 [dim s r2 dim s r2 min m] ms'; with -t: only lines with (min) r2 >= T, line + ' R<r2>'
 */
#include <pari/pari.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#undef ulong
#include "selmer_core.h"
#include "tables.h"

static u32 *spf; static u64 spf_lim;
#define NUMQ(t) (typ(t) == t_INT ? (t) : gel(t, 1))
#define DENQ(t) (typ(t) == t_INT ? gen_1 : gel(t, 2))
static seltab_t TT; static u32 qroff_h[QRMAX]; static u64 *qrbits_h;

static u64 rs_ = 0x9E3779B97F4A7C15ULL;
static double T_basis = 0, T_conic = 0, T_real = 0, T_fast = 0, T_padic = 0;
static inline double nowms(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return 1e3 * t.tv_sec + 1e-6 * t.tv_nsec; }
static inline u64 rnd64(void) { rs_ ^= rs_ << 13; rs_ ^= rs_ >> 7; rs_ ^= rs_ << 17; return rs_; }

/* ---- modular helpers (odd prime v < 2^62) ---- */
static u64 mm(u64 a, u64 b, u64 m) { return (u64)((unsigned __int128)a * b % m); }
static u64 pw(u64 a, u64 e, u64 m) { u64 r = 1 % m; a %= m; while (e) { if (e & 1) r = mm(r, a, m); a = mm(a, a, m); e >>= 1; } return r; }
static u64 inv(u64 a, u64 m) { return pw(a, m - 2, m); }
static u64 sqrtmod(u64 a, u64 p) {         /* Tonelli-Shanks, a a nonzero QR */
    if ((p & 3) == 3) return pw(a, (p + 1) / 4, p);
    u64 q = p - 1; int s = 0; while (!(q & 1)) { q >>= 1; s++; }
    u64 z = 2; while (legendre(z, p) != -1) z++;
    u64 c = pw(z, q, p), r = pw(a, (q + 1) / 2, p), t = pw(a, q, p); int m = s;
    while (t != 1) { int i = 0; u64 tt = t; while (tt != 1) { tt = mm(tt, tt, p); i++; } u64 b = c; for (int j = 0; j < m - i - 1; j++) b = mm(b, b, p); r = mm(r, b, p); c = mm(b, b, p); t = mm(t, c, p); m = i; }
    return r;
}

/* class of a nonzero element of Q_v: alpha = v(L) mod 2, u = unit info (odd v: 1 iff non-residue; v=2: unit mod 8; v=0: sign bit) */
typedef struct { int alpha, u; } lcls_t;

/* GEN of a product over a mask of the basis (bit 0 = -1, bit 1 = 2, bit q = B[q]) */
static GEN mask_gen(u64 m, const u64 *B, int n) {
    GEN r = gen_1;
    for (int q = 1; q < n; q++) if ((m >> q) & 1) r = mulii(r, utoi(B[q]));
    if (m & 1) r = negi(r);
    return r;
}

/* p-adic / real local point; returns 0 on success and fills cls[3] */
static int local_padic(GEN e[3], GEN l[3], GEN A[3][3], GEN sqv[3], ulong v, long prec, lcls_t cls[3], long maxtries) {
    pari_sp av0 = avma;
    GEN pv = utoipos(v);
    long K = prec / 2 + 4;
    GEN vK = powuu(v, K);
    pari_sp av = avma;                      /* reset point after pv, vK */
    for (long tries = 0; tries < maxtries; tries++, set_avma(av)) {
        GEN x; int m = (int)(rnd64() % 3);
        if (m == 0) { x = gdiv(randomi(vK), powuu(v, 2 * (rnd64() % 4))); }
        else {
            int c = (int)(rnd64() % 3); long j = (long)(rnd64() % (K + 4)) - 4;
            GEN off = addiu(randomi(vK), 1);
            x = gadd(e[c], j >= 0 ? mulii(powuu(v, j), off) : gdiv(off, powuu(v, -j)));
        }
        GEN u[3]; int ok = 1;
        for (int t = 0; t < 3 && ok; t++) {
            GEN w = gdiv(gsub(x, e[t]), l[t]);
            if (gequal0(w)) { ok = 0; break; }
            long val = gvaluation(w, pv); if (val & 1) { ok = 0; break; }
            GEN wp = gadd(w, zeropadic(pv, val + prec));
            if (!issquare(wp)) { ok = 0; break; }
            u[t] = gsqrt(wp, 0);
            if (rnd64() & 1) u[t] = gneg(u[t]);
        }
        if (!ok) continue;
        /* L_i = A0 sq + A1 u_j + A2 u_k with (j,k) = pairs */
        static const int pj[3] = {1, 0, 0}, pk[3] = {2, 2, 1};
        for (int i = 0; i < 3 && ok; i++) {
            GEN L = gadd(gadd(mulii(A[i][0], sqv[i]), gmul(A[i][1], u[pj[i]])), gmul(A[i][2], u[pk[i]]));
            if (gequal0(L)) { ok = 0; break; }
            long vl = valp(L); if (vl > prec / 2) { ok = 0; break; }
            GEN un = gmul(L, powis(pv, -vl));           /* unit part */
            cls[i].alpha = (int)(vl & 1);
            if (v == 2) { GEN t = padic_to_Q(un); ulong r = umodiu(NUMQ(t), 8); ulong d = umodiu(DENQ(t), 8); cls[i].u = (int)((r * d) & 7); /* d odd: d^-1 = d mod 8 */ }
            else { GEN t = padic_to_Q(un); ulong r = umodiu(NUMQ(t), v), d = umodiu(DENQ(t), v); cls[i].u = legendre(mm(r, inv(d, v), v), v) == -1; }
        }
        if (!ok) continue;
        set_avma(av0); return 0;
    }
    set_avma(av0); return -1;
}

static int local_real(GEN e[3], GEN l[3], GEN A[3][3], GEN sqv[3], lcls_t cls[3]) {
    pari_sp av = avma;
    long prec = nbits2prec(256);
    GEN es = sort(mkvec3(e[0], e[1], e[2]));
    GEN a = gel(es, 1), b = gel(es, 2), c = gel(es, 3);
    GEN cand[16]; int nc = 0;
    cand[nc++] = addii(addii(c, absi(c)), gen_1);
    cand[nc++] = subii(subii(a, absi(a)), gen_1);
    for (int k = 1; k <= 6; k++) { cand[nc++] = gadd(a, gdivgs(gmulgs(gsub(b, a), k), 7)); cand[nc++] = gadd(b, gdivgs(gmulgs(gsub(c, b), k), 7)); }
    for (int ci = 0; ci < nc; ci++) {
        GEN x = cand[ci], u[3]; int ok = 1;
        for (int t = 0; t < 3 && ok; t++) {
            GEN w = gdiv(gsub(x, e[t]), l[t]);
            if (gsigne(w) <= 0) { ok = 0; break; }
            u[t] = gsqrt(gtofp(w, prec), prec);
        }
        if (!ok) continue;
        static const int pj[3] = {1, 0, 0}, pk[3] = {2, 2, 1};
        for (int i = 0; i < 3 && ok; i++) {
            GEN t0 = gtofp(mulii(A[i][0], sqv[i]), prec), t1 = gmul(A[i][1], u[pj[i]]), t2 = gmul(A[i][2], u[pk[i]]);
            GEN L = gadd(gadd(t0, t1), t2);
            long em = maxss(maxss(expo(t0), expo(t1)), expo(t2));
            if (gequal0(L) || expo(L) < em - 150) { ok = 0; break; }
            cls[i].alpha = 0; cls[i].u = gsigne(L) < 0;
        }
        if (ok) { set_avma(av); return 0; }
    }
    set_avma(av); return -1;
}

/* fast odd-prime local point: v does not divide l1 l2 l3; returns 0 ok */
static int local_fast(u64 v, u64 er[3], u64 lr[3], u64 Ar[3][3], u64 sqr[3], lcls_t cls[3]) {
    if (v < 7) return -1;
    u64 li[3]; for (int t = 0; t < 3; t++) { if (lr[t] == 0) return -1; li[t] = inv(lr[t], v); }
    static const int pj[3] = {1, 0, 0}, pk[3] = {2, 2, 1};
    for (int tries = 0; tries < 256; tries++) {
        u64 x = rnd64() % v, u[3]; int ok = 1;
        for (int t = 0; t < 3 && ok; t++) {
            u64 w = mm((x + v - er[t]) % v, li[t], v);
            if (w == 0 || legendre(w, v) != 1) { ok = 0; break; }
            u[t] = sqrtmod(w, v);
        }
        if (!ok) continue;
        for (int i = 0; i < 3 && ok; i++) {
            u64 L = (mm(Ar[i][0], sqr[i], v) + mm(Ar[i][1], u[pj[i]], v) + mm(Ar[i][2], u[pk[i]], v)) % v;
            if (L == 0) { ok = 0; break; }
            cls[i].alpha = 0; cls[i].u = legendre(L, v) == -1;
        }
        if (ok) return 0;
    }
    return -1;
}

typedef struct { int dim, s, r2, err, sym; } casres_t;

static casres_t cassels_pair(int fam, int cidx, i64 p, i64 q, int *nfast, int *nslow) {
    casres_t R = {0, 0, 0, 0, 1};
    pari_sp av = avma;
    i64 F[MAXF]; int nf = family_forms(fam, p, q, F);
    if (nf <= 0) { R.err = 1; return R; }
    curve_t cv[3]; family_curves(fam, F, nf, cv);
    const curve_t *C = &cv[cidx];
    fact_t ff[MAXF];
    for (int m = 0; m < nf; m++) if (factor64(F[m] < 0 ? (u64)(-F[m]) : (u64)F[m], &ff[m], spf, spf_lim)) { R.err = 2; return R; }
    u64 B[MAXB]; int n; row_t P[2 * MAXB + 8]; int pv[2 * MAXB + 8], rank;
    double tb0 = nowms();
    int dim = sel2_dim_ex(C, ff, &TT, 0, B, &n, P, pv, &rank);
    T_basis += nowms() - tb0;
    if (dim < 0) { R.err = 3; return R; }
    R.dim = dim;
    /* RREF and kernel */
    for (int k = 0; k < rank; k++) for (int k2 = 0; k2 < rank; k2++) if (k2 != k) {
        int b = pv[k]; int has = b < 64 ? (int)((P[k2].lo >> b) & 1) : (int)((P[k2].hi >> (b - 64)) & 1);
        if (has) { P[k2].lo ^= P[k].lo; P[k2].hi ^= P[k].hi; }
    }
    int ispiv[128] = {0}; for (int k = 0; k < rank; k++) ispiv[pv[k]] = 1;
    u64 m1v[64], m2v[64]; int nb = 0;
    for (int f = 0; f < 2 * n; f++) if (!ispiv[f]) {
        int x[128] = {0}; x[f] = 1;
        for (int k = 0; k < rank; k++) { int bit = f < 64 ? (int)((P[k].lo >> f) & 1) : (int)((P[k].hi >> (f - 64)) & 1); x[pv[k]] = bit; }
        u64 a = 0, b = 0; for (int c = 0; c < n; c++) { if (x[c]) a |= 1ULL << c; if (x[n + c]) b |= 1ULL << c; }
        m1v[nb] = a; m2v[nb] = b; nb++;
    }
    if (nb != dim) { R.err = 4; return R; }
    /* root differences D1 = r, D2 = s, D3 = s - r as GEN; their class masks */
    GEN D[3]; u64 Dm[3];
    for (int t = 0; t < 3; t++) {
        GEN g = int2n(C->c2[t]);
        for (int m = 0; m < nf; m++) if (C->e[t][m]) g = mulii(g, powis(stoi(F[m]), C->e[t][m]));
        if (C->sgn[t] < 0) g = negi(g);
        if (signe(g) == 0) { R.err = 5; return R; }
        long v2 = C->c2[t];
        for (int m = 0; m < nf; m++) if (C->e[t][m]) for (int i = 0; i < ff[m].n; i++) if (ff[m].p[i] == 2) v2 += (long)C->e[t][m] * ff[m].e[i];
        u64 msk = (v2 & 1) ? 2 : 0;
        for (int qq = 2; qq < n; qq++) {
            long vv = 0;
            for (int m = 0; m < nf; m++) if (C->e[t][m]) for (int i = 0; i < ff[m].n; i++) if (ff[m].p[i] == B[qq]) vv += (long)C->e[t][m] * ff[m].e[i];
            if (vv & 1) msk |= 1ULL << qq;
        }
        if (signe(g) < 0) msk |= 1;
        D[t] = g; Dm[t] = msk;
    }
    GEN e[3] = { gen_0, D[0], D[1] };
    /* conic i: pair (j,k) = (2,3),(1,3),(1,2) (1-based); root difference e_k - e_j = D3, D2, D1 */
    static const int pj[3] = {1, 0, 0}, pk[3] = {2, 2, 1}, dt[3] = {2, 1, 0};
    GEN sq[3]; GEN c0[3];
    for (int i = 0; i < 3; i++) { c0[i] = mask_gen(Dm[dt[i]], B, n); GEN qq = diviiexact(D[dt[i]], c0[i]); if (signe(qq) < 0) { R.err = 6; return R; } sq[i] = sqrti(qq); if (!equalii(sqri(sq[i]), qq)) { R.err = 6; return R; } }
    /* drop the image of E(Q)[2] (the pairing vanishes on it): keep a complement of span{T1, T2} in the kernel basis.
       Kummer images: T1 = (e1,0) -> (D1 D2, -D1), T2 = (e2,0) -> (D1, -D1 D3). */
    {
        row_t piv[128]; u64 pm[2] = {0, 0}; int rk = 0;
        u64 tv1[2] = { Dm[0] ^ Dm[1], Dm[0] ^ 1 }, tv2[2] = { Dm[0], Dm[0] ^ Dm[2] ^ 1 };
        u64 *tvs[2] = { tv1, tv2 };
        for (int k = 0; k < 2; k++) { row_t y; y.lo = tvs[k][0] | (tvs[k][1] << n); y.hi = tvs[k][1] >> (64 - n); insert_row_t(piv, pm, &rk, y, (row_t *)0, (int *)0); }
        int nkeep = 0;
        for (int a = 0; a < nb; a++) {
            int before = rk; row_t y; y.lo = m1v[a] | (m2v[a] << n); y.hi = m2v[a] >> (64 - n);
            insert_row_t(piv, pm, &rk, y, (row_t *)0, (int *)0);
            if (rk > before) { m1v[nkeep] = m1v[a]; m2v[nkeep] = m2v[a]; nkeep++; }
        }
        nb = nkeep;
    }
    /* precompute Legendre masks of the basis at each odd place */
    u64 nonresm[MAXB]; for (int j = 2; j < n; j++) { u64 m = 0; u64 v = B[j]; if ((v & 3) == 3) m |= 1; { int r8 = (int)(v & 7); if (r8 == 3 || r8 == 5) m |= 2; } for (int qq = 2; qq < n; qq++) if (qq != j && legendre(B[qq], v) == -1) m |= 1ULL << qq; nonresm[j] = m; }
    u64 H[64][3];
    long prec_v[MAXB], prec_lo[MAXB];
    for (int j = 1; j < n; j++) { long vv = 0; for (int t = 0; t < 3; t++) { long x = gvaluation(D[t], utoipos(B[j])); vv += x; } prec_v[j] = 2 * vv + 40; prec_lo[j] = vv + 12; }
    for (int a = 0; a < nb; a++) {
        pari_sp av2 = avma;
        u64 lm[3] = { m1v[a], m2v[a], m1v[a] ^ m2v[a] };
        GEN l[3]; for (int t = 0; t < 3; t++) l[t] = mask_gen(lm[t], B, n);
        GEN A[3][3];
        double tc0 = nowms();
        for (int i = 0; i < 3; i++) {
            int j = pj[i], k = pk[i];
            GEN G = mkmat3(mkcol3(negi(c0[i]), gen_0, gen_0), mkcol3(gen_0, l[j], gen_0), mkcol3(gen_0, gen_0, negi(l[k])));
            u64 pm = Dm[dt[i]] | lm[j] | lm[k]; GEN Pl = cgetg(1, t_VEC);
            for (int qq = 1; qq < n; qq++) if ((pm >> qq) & 1) Pl = vec_append(Pl, utoipos(B[qq]));
            GEN sol = (lg(Pl) > 1) ? qfsolve(mkvec2(G, Pl)) : qfsolve(G);
            if (typ(sol) != t_COL) { R.err = 7; set_avma(av); return R; }
            sol = Q_primpart(sol);
            A[i][0] = mulii(negi(c0[i]), gel(sol, 1)); A[i][1] = mulii(l[j], gel(sol, 2)); A[i][2] = mulii(negi(l[k]), gel(sol, 3));
        }
        T_conic += nowms() - tc0;
        H[a][0] = H[a][1] = H[a][2] = 0;
        /* infinity */
        lcls_t cls[3];
        double tr0 = nowms();
        if (local_real(e, l, A, sq, cls)) { R.err = 8; set_avma(av); return R; }
        T_real += nowms() - tr0;
        for (int i = 0; i < 3; i++) if (cls[i].u) H[a][i] ^= 1;   /* (L,-1)_inf */
        for (int j = 1; j < n; j++) {
            u64 v = B[j]; int done = 0;
            if (v > 2 && !(((lm[0] | lm[1] | lm[2]) >> j) & 1)) {
                u64 er[3] = { 0, umodiu(e[1], v), umodiu(e[2], v) }, lr[3], Ar[3][3], sqr[3];
                for (int t = 0; t < 3; t++) lr[t] = umodiu(l[t], v);
                for (int i = 0; i < 3; i++) { for (int c = 0; c < 3; c++) Ar[i][c] = umodiu(A[i][c], v); sqr[i] = umodiu(sq[i], v); }
                double tf0 = nowms();
                if (!local_fast(v, er, lr, Ar, sqr, cls)) { done = 1; (*nfast)++; }
                T_fast += nowms() - tf0;
            }
            if (!done) {
                double tp0 = nowms();
                /* low precision first (enough in almost all cases), then the safe high precision */
                if (local_padic(e, l, A, sq, (ulong)v, prec_lo[j], cls, 3000) &&
                    local_padic(e, l, A, sq, (ulong)v, prec_v[j], cls, 400000)) { R.err = 9; set_avma(av); return R; }
                T_padic += nowms() - tp0;
                (*nslow)++;
            }
            for (int i = 0; i < 3; i++) {
                int al = cls[i].alpha; u64 hm = 0;
                if (v == 2) {
                    int u = cls[i].u, e3 = (u & 3) == 3, w3 = (u == 3 || u == 5);
                    if (e3) hm |= 1;                    /* (L,-1)_2 = (-1)^eps(u) */
                    if (w3) hm |= 2;                    /* (L,2)_2 = (-1)^omega(u) */
                    for (int qq = 2; qq < n; qq++) { u64 b = B[qq]; int be3 = (b & 3) == 3, bw = ((b & 7) == 3 || (b & 7) == 5); if ((e3 & be3) ^ (al & bw)) hm |= 1ULL << qq; }
                } else {
                    if (al) hm = nonresm[j];            /* (L,b)_v = (b|v)^alpha for b != v, incl. -1 and 2 */
                    int self = (al && (v & 3) == 3) ^ cls[i].u;   /* (L,v)_v = (-1)^{alpha eps(v)} (u|v) */
                    hm = (hm & ~(1ULL << j)) | ((u64)self << j);
                }
                H[a][i] ^= hm;
            }
        }
        set_avma(av2);
    }
    /* pairing matrix and rank */
    u64 Crow[64];
    for (int a = 0; a < nb; a++) {
        u64 r = 0;
        for (int b = 0; b < nb; b++) {
            u64 lb[3] = { m1v[b], m2v[b], m1v[b] ^ m2v[b] };
            int c = par64(H[a][0] & lb[0]) ^ par64(H[a][1] & lb[1]) ^ par64(H[a][2] & lb[2]);
            if (c) r |= 1ULL << b;
        }
        Crow[a] = r;
    }
    for (int a = 0; a < nb; a++) { if ((Crow[a] >> a) & 1) R.sym = 0; for (int b = 0; b < nb; b++) if (((Crow[a] >> b) & 1) != ((Crow[b] >> a) & 1)) R.sym = 0; }
    int rk = 0;
    for (int col = 0; col < nb; col++) {
        int piv = -1; for (int r = rk; r < nb; r++) if ((Crow[r] >> col) & 1) { piv = r; break; }
        if (piv < 0) continue;
        u64 t = Crow[piv]; Crow[piv] = Crow[rk]; Crow[rk] = t;
        for (int r = 0; r < nb; r++) if (r != rk && ((Crow[r] >> col) & 1)) Crow[r] ^= Crow[rk];
        rk++;
    }
    R.s = rk; R.r2 = dim - 2 - rk;
    set_avma(av);
    return R;
}

int main(int argc, char **argv) {
    int target = -1000;
    for (int i = 1; i < argc; i++) if (!strcmp(argv[i], "-t") && i + 1 < argc) target = atoi(argv[++i]);
    pari_init(1L << 27, 0);
    setrand(utoipos(12345));
    spf_lim = 1u << 24; spf = (u32 *)calloc(spf_lim + 1, 4);
    for (u64 i = 2; i <= spf_lim; i++) if (!spf[i]) { spf[i] = (u32)i; if (i * i <= spf_lim) for (u64 j = i * i; j <= spf_lim; j += i) if (!spf[j]) spf[j] = (u32)i; }
    qrbits_h = (u64 *)calloc(QR_WORDS, 8); qr_build(qroff_h, qrbits_h, QR_WORDS);
    TT.tab2 = &tab2[0][0]; TT.tab3 = &tab3[0][0]; TT.qroff = qroff_h; TT.qrbits = qrbits_h; TT.tabw2 = 0;
    char line[4096], fs[32]; long long p, q; long nin = 0, nout = 0, nfast = 0, nslow = 0;
    double ttot = 0;
    while (fgets(line, sizeof line, stdin)) {
        if (sscanf(line, "%31s %lld %lld", fs, &p, &q) != 3) continue;
        int fam = !strcmp(fs, "z8") ? FAM_Z8 : !strcmp(fs, "z2z6") ? FAM_Z2Z6 : !strcmp(fs, "z12") ? FAM_Z12 : !strcmp(fs, "z2z8") ? FAM_Z2Z8 : -1;
        if (fam < 0) continue;
        nin++;
        int ncur = (fam == FAM_Z2Z8) ? 3 : 1;      /* z8: E1 only; z2z8: the three full-2-torsion curves (min r2) */
        int nf_ = 0, ns_ = 0;
        struct timespec t0, t1; clock_gettime(CLOCK_MONOTONIC, &t0);
        casres_t Rs[3]; int minr2 = 1000, err = 0, asym = 0, done_c = 0;
        for (int c = 0; c < ncur; c++) {
            Rs[c] = cassels_pair(fam, c, p, q, &nf_, &ns_); done_c++;
            if (Rs[c].err) { err = Rs[c].err; break; }
            if (!Rs[c].sym) asym = 1;
            if (Rs[c].r2 < minr2) minr2 = Rs[c].r2;
            if (target != -1000 && minr2 < target) break;      /* min already below target */
        }
        clock_gettime(CLOCK_MONOTONIC, &t1);
        double ms = 1e3 * ((t1.tv_sec - t0.tv_sec) + 1e-9 * (t1.tv_nsec - t0.tv_nsec)); ttot += ms;
        nfast += nf_; nslow += ns_;
        size_t L = strlen(line); while (L && (line[L - 1] == '\n' || line[L - 1] == '\r')) line[--L] = 0;
        if (target == -1000) {
            printf("%s %lld %lld", fs, p, q);
            for (int c = 0; c < done_c; c++) printf(" %d %d %d", Rs[c].dim, Rs[c].s, Rs[c].r2);
            if (ncur > 1) printf(" min %d", err ? 99 : minr2);
            printf(" %.3f%s%s\n", ms, err ? " ERR" : "", asym ? " ASYM" : "");
        }
        else if (err || minr2 >= target) { printf("%s R%d%s\n", line, err ? 99 : minr2, err ? " ERR" : ""); nout++; }
        if (err) fprintf(stderr, "error %d on %s %lld %lld\n", err, fs, p, q);
    }
    fprintf(stderr, "casc: %ld pairs, %.3f ms/pair, local points fast %ld / padic %ld%s\n", nin, nin ? ttot / nin : 0, nfast, nslow, target == -1000 ? "" : "");
    if (target != -1000) fprintf(stderr, "casc: kept %ld\n", nout);
    if (nin) fprintf(stderr, "casc: ms/pair: basis %.3f conics %.3f real %.3f fast %.3f padic %.3f\n", T_basis / nin, T_conic / nin, T_real / nin, T_fast / nin, T_padic / nin);
    pari_close();
    return 0;
}
