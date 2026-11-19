/* selmer_core.h -- exact dim Sel2 for the full-2-torsion curves of the Kubert families.
 * Shared, bit for bit, by the CPU build (selmer.c, gcc) and the GPU build (selmer_gpu.cu, nvcc).
 *
 * Curve y^2 = x (x - r)(x - s), roots e = (0, r, s); D1 = r = e2-e1, D2 = s = e3-e1, D3 = s-r = e3-e2,
 * each D_t = sgn * 2^c2 * prod_m F_m^e[t][m]  with F_m the cusp forms (int64 values).
 * Sel2 = kernel of the F2 matrix of local conditions  P_v . res_v(d1,d2) = 0,  (d1,d2) = classes of (x-e1, x-e2),
 * P_v a basis of the annihilator (dot product) of the local image W_v:
 *   v = 2: table tab2 (key: min valuation mod 2, deeper difference, gap capped, units mod 16)  [oracle-built]
 *   v = 3: table tab3 (same key, units mod 9)                                                   [oracle-built]
 *   v >= 5: rules (split / non-split by the Legendre symbol of e_i-e_k, parity of v_l(e_i-e_j))
 *   v = inf: the class of the egg component.
 */
#pragma once
#define SEL_SINGULAR 255
#define SEL_FACTFAIL 254
#define SEL_ERROR 253
#include <stdint.h>
#include <math.h>

#ifdef __CUDACC__
#define HD __host__ __device__ __forceinline__
#define HDN __host__ __device__ __noinline__
#else
#define HD static inline
#define HDN static
#endif

typedef uint64_t u64;
typedef int64_t i64;
typedef uint32_t u32;

/* ---------------- 64-bit Montgomery arithmetic (n odd, n < 2^62) ---------------- */
HD u64 mulhi64(u64 a, u64 b) {
#ifdef __CUDA_ARCH__
    return __umul64hi(a, b);
#else
    return (u64)(((unsigned __int128)a * b) >> 64);
#endif
}
/* 64-bit % and / with a 32-bit fast path (the GPU emulates 64-bit division in ~100 instructions) */
HD u64 umod(u64 a, u64 m) { return ((a | m) >> 32) ? a % m : (u64)((u32)a % (u32)m); }
HD u64 udiv(u64 a, u64 m) { return ((a | m) >> 32) ? a / m : (u64)((u32)a / (u32)m); }
typedef struct { u64 n, ninv, r2, one; } mont_t;     /* ninv = -n^{-1} mod 2^64, one = 2^64 mod n, r2 = 2^128 mod n */
HD u64 mont_redc(u64 hi, u64 lo, const mont_t *M) {
    u64 m = lo * M->ninv;
    u64 t = hi + mulhi64(m, M->n) + (lo != 0);
    return t >= M->n ? t - M->n : t;
}
HD u64 mont_mul(u64 a, u64 b, const mont_t *M) { return mont_redc(mulhi64(a, b), a * b, M); }
HD void mont_init(mont_t *M, u64 n) {
    M->n = n;
    u64 inv = n;                                   /* Newton: inv = n^{-1} mod 2^64 */
    for (int i = 0; i < 5; i++) inv *= 2 - n * inv;
    M->ninv = (u64)0 - inv;
    u64 r = ((u64)0 - n) % n;                      /* 2^64 mod n */
    M->one = r;
    /* r2 = 2^128 mod n by doubling r 64 times */
    u64 x = r;
    for (int i = 0; i < 64; i++) { x = (x >= n - x) ? x - (n - x) : x + x; }
    M->r2 = x;
}
HD u64 mont_to(u64 a, const mont_t *M) { return mont_mul(a % M->n, M->r2, M); }
HD u64 mont_from(u64 a, const mont_t *M) { return mont_redc(0, a, M); }

/* plain modular multiply for small moduli (< 2^32) and generic via Montgomery */
HD u64 mulmod_any(u64 a, u64 b, u64 m) {
    if (m < (1ULL << 16)) return (u64)(((u32)umod(a, m) * (u32)umod(b, m)) % (u32)m);
    if (m < (1ULL << 32)) return umod(a, m) * umod(b, m) % m;
#ifdef __CUDA_ARCH__
    /* 128-bit remainder without __int128: binary method (rare path) */
    u64 r = 0; a %= m; b %= m;
    while (b) { if (b & 1) { r = (r >= m - a) ? r - (m - a) : r + a; } a = (a >= m - a) ? a - (m - a) : a + a; b >>= 1; }
    return r;
#else
    return (u64)((unsigned __int128)a * b % m);
#endif
}
HD u64 powmod_any(u64 a, u64 e, u64 m) {
    u64 r = 1 % m; a %= m;
    while (e) { if (e & 1) r = mulmod_any(r, a, m); a = mulmod_any(a, a, m); e >>= 1; }
    return r;
}

/* deterministic Miller-Rabin for n < 2^64 (n odd > 37) */
HD int mr_prime(u64 n) {
    mont_t M; mont_init(&M, n);
    u64 d = n - 1; int s = 0; while (!(d & 1)) { d >>= 1; s++; }
    const u64 bases[7] = {2, 325, 9375, 28178, 450775, 9780504, 1795265022};
    u64 mone = M.n - M.one;  /* -1 in Montgomery form */
    for (int i = 0; i < 7; i++) {
        u64 a = bases[i] % n; if (a == 0) continue;
        u64 x = M.one, b = mont_to(a, &M), e = d;
        while (e) { if (e & 1) x = mont_mul(x, b, &M); b = mont_mul(b, b, &M); e >>= 1; }
        if (x == M.one || x == mone) continue;
        int ok = 0;
        for (int r = 1; r < s; r++) { x = mont_mul(x, x, &M); if (x == mone) { ok = 1; break; } }
        if (!ok) return 0;
    }
    return 1;
}
HD int is_prime64(u64 n) {
    if (n < 2) return 0;
    const u32 sp[12] = {2,3,5,7,11,13,17,19,23,29,31,37};
    for (int i = 0; i < 12; i++) if (n % sp[i] == 0) return n == sp[i];
    if (n < 1369) return 1;
    return mr_prime(n);
}
HD u64 gcd64(u64 a, u64 b) { while (b) { u64 t = a % b; a = b; b = t; } return a; }
/* Pollard-Brent rho, Montgomery form; returns a nontrivial factor or 0 after the iteration cap */
HD u64 rho64(u64 n, u64 cap) {
    if (!(n & 1)) return 2;
    mont_t M; mont_init(&M, n);
    for (u64 c0 = 1; c0 < 20; c0++) {
        u64 c = mont_to(c0, &M);
        u64 y = mont_to(2, &M), x = y, ys = y, q = M.one, g = 1, r = 1, it = 0;
        const u64 m = 64;
        do {
            x = y;
            for (u64 i = 0; i < r; i++) { y = mont_mul(y, y, &M) + c; if (y >= n) y -= n; }
            u64 k = 0;
            do {
                ys = y;
                u64 lim = (m < r - k) ? m : r - k;
                for (u64 i = 0; i < lim; i++) { y = mont_mul(y, y, &M) + c; if (y >= n) y -= n; q = mont_mul(q, x > y ? x - y : y - x, &M); }
                g = gcd64(q, n); k += m; it += lim;
            } while (k < r && g == 1);
            r <<= 1;
        } while (g == 1 && it < cap);
        if (g == n) {
            g = 1;
            for (u64 i = 0; i < 4 * r && g == 1; i++) { ys = mont_mul(ys, ys, &M) + c; if (ys >= n) ys -= n; g = gcd64(x > ys ? x - ys : ys - x, n); }
        }
        if (g != 1 && g != n) return g;
        if (it >= cap) return 0;
    }
    return 0;
}

HD int ctz64(u64 x) {
#ifdef __CUDA_ARCH__
    return __ffsll((long long)x) - 1;
#else
    return __builtin_ctzll(x);
#endif
}
HD int ctz32(u32 x) {
#ifdef __CUDA_ARCH__
    return __ffs((int)x) - 1;
#else
    return __builtin_ctz(x);
#endif
}
/* Jacobi symbol (a/n), n odd: binary algorithm, no division */
HD int jacobi64(u64 a, u64 n) {
    int t = 1;
    if (a >= n) a %= n;
    while (a) {
        int z = ctz64(a);
        a >>= z;
        if ((z & 1) && ((n & 7) == 3 || (n & 7) == 5)) t = -t;
        if (a < n) { u64 tmp = a; a = n; n = tmp; if ((a & 3) == 3 && (n & 3) == 3) t = -t; }
        a -= n;
    }
    return n == 1 ? t : 0;
}
HD int jacobi32(u32 a, u32 n) {
    int t = 1;
    if (a >= n) a %= n;
    while (a) {
        int z = ctz32(a);
        a >>= z;
        if ((z & 1) && ((n & 7) == 3 || (n & 7) == 5)) t = -t;
        if (a < n) { u32 tmp = a; a = n; n = tmp; if ((a & 3) == 3 && (n & 3) == 3) t = -t; }
        a -= n;
    }
    return n == 1 ? t : 0;
}
HD int legendre(u64 a, u64 l) { return (l < (1ULL << 32)) ? jacobi32((u32)umod(a, l), (u32)l) : jacobi64(a, l); }


/* ---------------- factorisation ---------------- */
#define MAXPF 16
typedef struct { u64 p[MAXPF]; uint8_t e[MAXPF]; int n; } fact_t;

HD int fact_add(fact_t *f, u64 p, int e) {
    for (int i = 0; i < f->n; i++) if (f->p[i] == p) { f->e[i] += e; return 0; }
    if (f->n >= MAXPF) return -1;
    f->p[f->n] = p; f->e[f->n] = (uint8_t)e; f->n++; return 0;
}
/* factor n >= 1 using the spf table (n <= lim), trial division and rho; returns 0 ok, <0 failure */
HDN int factor64(u64 n, fact_t *f, const u32 *spf, u64 lim) {
    f->n = 0;
    if (n <= lim) {
        while (n > 1) { u64 p = spf[n]; int e = 0; while (umod(n, p) == 0) { n = udiv(n, p); e++; } if (fact_add(f, p, e)) return -1; }
        return 0;
    }
    /* trial division by small primes */
    const u32 tp[25] = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97};
    for (int i = 0; i < 25; i++) {
        u32 p = tp[i];
        if (umod(n, p) == 0) { int e = 0; while (umod(n, p) == 0) { n = udiv(n, p); e++; } if (fact_add(f, p, e)) return -1; }
    }
    /* stack of cofactors */
    u64 st[16]; int sp = 0;
    if (n > 1) st[sp++] = n;
    while (sp) {
        u64 m = st[--sp];
        if (m <= lim) {
            while (m > 1) { u64 p = spf[m]; int e = 0; while (umod(m, p) == 0) { m = udiv(m, p); e++; } if (fact_add(f, p, e)) return -1; }
            continue;
        }
        if (m < 9409 || is_prime64(m)) { if (fact_add(f, m, 1)) return -1; continue; }   /* no factor < 97 */
        u64 d = rho64(m, 1ULL << 22);
        if (!d) return -2;
        if (sp + 2 > 16) return -3;
        st[sp++] = d; st[sp++] = m / d;
    }
    return 0;
}

/* ---------------- curve description ---------------- */
#define MAXF 8
typedef struct {
    int nf;
    i64 F[MAXF];
    int sgn[3];
    int c2[3];
    int8_t e[3][MAXF];
} curve_t;

enum { FAM_Z8 = 0, FAM_Z2Z6 = 1, FAM_Z12 = 2, FAM_Z2Z8 = 3, FAM_Z10 = 4, FAM_Z9 = 5 };

/* forms of a family (shared by its curves); returns nf, or 0 if a form vanishes (singular fibre) */
HD int family_forms(int fam, i64 p, i64 q, i64 *F) {
    int nf = 0;
    switch (fam) {
    case FAM_Z8:   /* p, p-q, u=2p-q, q, Q=8p^2-8pq+q^2 */
        F[0] = p; F[1] = p - q; F[2] = 2 * p - q; F[3] = q; F[4] = 8 * p * p - 8 * p * q + q * q; nf = 5; break;
    case FAM_Z2Z6: /* q, p-3q, p-5q, p+3q, p-9q, p-q */
        F[0] = q; F[1] = p - 3 * q; F[2] = p - 5 * q; F[3] = p + 3 * q; F[4] = p - 9 * q; F[5] = p - q; nf = 6; break;
    case FAM_Z12:  /* q, u=2p-q, A2=2p^2-2pq+q^2, A6=6p^2-6pq+q^2, p, p-q, A3=3p^2-3pq+q^2 */
        F[0] = q; F[1] = 2 * p - q; F[2] = 2 * p * p - 2 * p * q + q * q; F[3] = 6 * p * p - 6 * p * q + q * q;
        F[4] = p; F[5] = p - q; F[6] = 3 * p * p - 3 * p * q + q * q; nf = 7; break;
    case FAM_Z2Z8: /* p, q, 2p+q, 4p+q, A=8p^2+4pq+q^2, B=8p^2-q^2, C=8p^2+8pq+q^2 */
        F[0] = p; F[1] = q; F[2] = 2 * p + q; F[3] = 4 * p + q; F[4] = 8 * p * p + 4 * p * q + q * q;
        F[5] = 8 * p * p - q * q; F[6] = 8 * p * p + 8 * p * q + q * q; nf = 7; break;
    case FAM_Z10:  /* p, p-q, P2=p^2-3pq+q^2, q, u=2p-q, F4=4p^2-2pq-q^2  (B = 16 p^5 (p-q)^5 P2, B' = q^5 u^5 F4) */
        F[0] = p; F[1] = p - q; F[2] = p * p - 3 * p * q + q * q; F[3] = q; F[4] = 2 * p - q; F[5] = 4 * p * p - 2 * p * q - q * q; nf = 6; break;
    case FAM_Z9:   /* p, q, p-q (width 9), A9 = p^2-pq+q^2 (width 3), C3 = p^3-6p^2q+3pq^2+q^3 (width 1) */
        F[0] = p; F[1] = q; F[2] = p - q; F[3] = p * p - p * q + q * q; F[4] = p * p * p - 6 * p * p * q + 3 * p * q * q + q * q * q; nf = 5; break;
    default: return -1;
    }
    for (int i = 0; i < nf; i++) if (F[i] == 0) return 0;
    return nf;
}

HD void curve_set(curve_t *c, int t, int sgn, int c2) { c->sgn[t] = sgn; c->c2[t] = c2; }
/* the full-2-torsion curves of the family; returns their number */
HD int family_curves(int fam, const i64 *F, int nf, curve_t *cv) {
    for (int c = 0; c < 3; c++) {
        cv[c].nf = nf;
        for (int m = 0; m < MAXF; m++) { cv[c].F[m] = m < nf ? F[m] : 0; for (int t = 0; t < 3; t++) cv[c].e[t][m] = 0; }
        for (int t = 0; t < 3; t++) { cv[c].sgn[t] = 1; cv[c].c2[t] = 0; }
    }
    switch (fam) {
    case FAM_Z8:
        /* E1: D1 = u^4, D2 = q^2 Q, D3 = -16 p^2 (p-q)^2 */
        cv[0].e[0][2] = 4; cv[0].e[1][3] = 2; cv[0].e[1][4] = 1;
        curve_set(&cv[0], 2, -1, 4); cv[0].e[2][0] = 2; cv[0].e[2][1] = 2;
        /* E2: D1 = q^4, D2 = Q^2, D3 = 16 p (p-q) u^2 */
        cv[1].e[0][3] = 4; cv[1].e[1][4] = 2;
        curve_set(&cv[1], 2, 1, 4); cv[1].e[2][0] = 1; cv[1].e[2][1] = 1; cv[1].e[2][2] = 2;
        return 2;
    case FAM_Z2Z6:
        /* D1 = -128 q^3 (p-3q), D2 = -(p-5q)^3 (p+3q), D3 = -(p-9q)(p-q)^3 */
        curve_set(&cv[0], 0, -1, 7); cv[0].e[0][0] = 3; cv[0].e[0][1] = 1;
        curve_set(&cv[0], 1, -1, 0); cv[0].e[1][2] = 3; cv[0].e[1][3] = 1;
        curve_set(&cv[0], 2, -1, 0); cv[0].e[2][4] = 1; cv[0].e[2][5] = 3;
        return 1;
    case FAM_Z12:
        /* D1 = q^2 u^6, D2 = A2^3 A6, D3 = 16 p^3 (p-q)^3 A3 */
        cv[0].e[0][0] = 2; cv[0].e[0][1] = 6; cv[0].e[1][2] = 3; cv[0].e[1][3] = 1;
        curve_set(&cv[0], 2, 1, 4); cv[0].e[2][4] = 3; cv[0].e[2][5] = 3; cv[0].e[2][6] = 1;
        return 1;
    case FAM_Z2Z8:
        /* #1: D1 = 256 p^4 (2p+q)^4, D2 = B A^2 C, D3 = -q^4 (4p+q)^4 */
        curve_set(&cv[0], 0, 1, 8); cv[0].e[0][0] = 4; cv[0].e[0][2] = 4;
        cv[0].e[1][5] = 1; cv[0].e[1][4] = 2; cv[0].e[1][6] = 1;
        curve_set(&cv[0], 2, -1, 0); cv[0].e[2][1] = 4; cv[0].e[2][3] = 4;
        /* #2: D1 = A^4, D2 = B^2 C^2, D3 = -64 q^2 p^2 (2p+q)^2 (4p+q)^2 */
        cv[1].e[0][4] = 4; cv[1].e[1][5] = 2; cv[1].e[1][6] = 2;
        curve_set(&cv[1], 2, -1, 6); cv[1].e[2][1] = 2; cv[1].e[2][0] = 2; cv[1].e[2][2] = 2; cv[1].e[2][3] = 2;
        /* #3: D1 = C^4, D2 = B^4, D3 = -32 q p (2p+q)(4p+q) A^2 */
        cv[2].e[0][6] = 4; cv[2].e[1][5] = 4;
        curve_set(&cv[2], 2, -1, 5); cv[2].e[2][1] = 1; cv[2].e[2][0] = 1; cv[2].e[2][2] = 1; cv[2].e[2][3] = 1; cv[2].e[2][4] = 2;
        return 3;
    }
    return 0;
}

/* ---------------- Selmer rank ---------------- */
#define MAXB 56
#define T2_K 16
#define T2_CAP 8
#define T3_K 9
#define T3_CAP 4
typedef struct { u64 lo, hi; } row_t;
HD void row_set(row_t *r, int c) { if (c < 64) r->lo |= 1ULL << c; else r->hi |= 1ULL << (c - 64); }
HD int par64(u64 x) {
#ifdef __CUDA_ARCH__
    return __popcll(x) & 1;
#else
    return __builtin_parityll(x);
#endif
}

/* unit part of D_t modulo M at prime l (M = l, 9 or 16), sign included */
HD u64 unit_res(const curve_t *C, int t, u64 l, u64 M, const int *vf) {
    u64 r = C->sgn[t] > 0 ? 1 : M - 1;
    if (l != 2) r = mulmod_any(r, powmod_any(2, (u64)C->c2[t], M), M);
    for (int m = 0; m < C->nf; m++) {
        int e = C->e[t][m]; if (!e) continue;
        i64 F = C->F[m]; int neg = F < 0; u64 a = neg ? (u64)(-F) : (u64)F;
        for (int k = 0; k < vf[m]; k++) a /= l;
        u64 am = a % M; if (neg) am = (M - am) % M;
        r = mulmod_any(r, powmod_any(am, (u64)e, M), M);
    }
    return r;
}

HD void annih4(int g1, int g2, int out[2]) {
    int k = 0;
    for (int a = 1; a < 16 && k < 2; a++) {
        if (par64((u64)(a & g1)) || par64((u64)(a & g2))) continue;
        if (k == 1 && a == out[0]) continue;
        out[k++] = a;
    }
}

/* dim Sel2 of curve C whose forms are factored in ff[]; returns dim, or a negative error code */
/* tables handed to the core (host or device pointers) */
#define QRMAX 4096
typedef struct {
    const uint8_t *tab2, *tab3;
    const u32 *qroff;     /* qroff[b] = bit offset of the QR bitmap of the odd prime b < QRMAX */
    const u64 *qrbits;    /* bit r of bitmap(b) set iff r is a nonzero square mod b */
    const uint8_t *tabw2; /* z8 parity cut: local root number at 2 by the tab2 key (tabw2.h); may be NULL */
} seltab_t;

/* 1 iff b is a non-residue mod l (b, l distinct odd primes) */
HD int nonres(u64 b, u64 l, const seltab_t *T) {
    if (l < QRMAX) { u64 r = umod(b, l), o = T->qroff[l] + r; return !((T->qrbits[o >> 6] >> (o & 63)) & 1); }
    if (b < QRMAX) {
        u64 r = umod(l, b), o = T->qroff[b] + r;
        int nr = !((T->qrbits[o >> 6] >> (o & 63)) & 1);
        return nr ^ (int)(((b & 3) == 3) && ((l & 3) == 3));
    }
    return legendre(b, l) == -1;
}

HD u64 mulm(u64 a, u64 b, u64 m) { return mulmod_any(a, b, m); }

/* dim Sel2 of curve C whose forms are factored in ff[]; returns dim, or a negative error code */
/* echelon insertion with a pivot table: piv[c] = row whose lowest set bit is column c (pmask marks used columns).
   Each reduction step clears the lowest set bit of y, so the loop runs at most popcount-many steps. */
HD void insert_row_t(row_t *piv, u64 *pmask, int *rank, row_t y, row_t *P, int *pv) {
    while (y.lo | y.hi) {
        int c = y.lo ? ctz64(y.lo) : 64 + ctz64(y.hi);
        int used = c < 64 ? (int)((pmask[0] >> c) & 1) : (int)((pmask[1] >> (c - 64)) & 1);
        if (!used) {
            piv[c] = y; if (c < 64) pmask[0] |= 1ULL << c; else pmask[1] |= 1ULL << (c - 64);
            if (P) { P[*rank] = y; pv[*rank] = c; }
            (*rank)++; return;
        }
        y.lo ^= piv[c].lo; y.hi ^= piv[c].hi;
    }
}
HD void insert_row(row_t *P, int *pv, int *rank, row_t y) {
    for (int k = 0; k < *rank; k++) {
        int b = pv[k];
        if (b < 64 ? ((y.lo >> b) & 1) : ((y.hi >> (b - 64)) & 1)) { y.lo ^= P[k].lo; y.hi ^= P[k].hi; }
    }
    if (y.lo || y.hi) { pv[*rank] = y.lo ? ctz64(y.lo) : 64 + ctz64(y.hi); P[*rank] = y; (*rank)++; }
}
/* optional outputs (host use, pass NULL on the device): Bout[0..n) basis of Q(S,2) (-1, 2, odd primes), *nout = n,
   Pout/pvout/rankout = echelon rows (2n columns: d1 bits 0..n-1, d2 bits n..2n-1) and their pivot columns. */
HDN int sel2_dim_ex(const curve_t *C, const fact_t *ff, const seltab_t *T, int stop_dim, u64 *Bout, int *nout, row_t *Pout, int *pvout, int *rankout) {
    u64 B[MAXB]; int n = 2;
    B[0] = (u64)-1; B[1] = 2;
    int nf = C->nf;
    int used[MAXF];
    for (int m = 0; m < nf; m++) {
        used[m] = (C->e[0][m] | C->e[1][m] | C->e[2][m]) != 0;
        if (!used[m]) continue;
        for (int i = 0; i < ff[m].n; i++) {
            u64 p = ff[m].p[i]; if (p == 2) continue;
            int seen = 0; for (int j = 2; j < n; j++) if (B[j] == p) { seen = 1; break; }
            if (!seen) { if (n >= MAXB) return -10; B[n++] = p; }
        }
    }
#if defined(STAGE) && STAGE == 1
    return n;
#endif
    row_t R0 = {0, 0};
    u32 mk[MAXB]; uint8_t meta[MAXB];
    int sg[3];
    for (int t = 0; t < 3; t++) {
        int s = C->sgn[t];
        for (int m = 0; m < nf; m++) if ((C->e[t][m] & 1) && C->F[m] < 0) s = -s;
        sg[t] = s;
    }
    {   /* real place */
        int zero_min = (sg[0] > 0 && sg[1] > 0), r_min = (sg[0] < 0 && sg[2] > 0);
        if (!r_min) row_set(&R0, 0);
        if (!zero_min) row_set(&R0, n);
    }
    /* 2-adic class masks of the basis: bit q of m1 = (-1)-component, m5 = 5-component */
    u64 m1 = 1, m5 = 0;
    for (int q = 2; q < n; q++) { int r8 = (int)(B[q] & 7); if ((r8 & 3) == 3) m1 |= 1ULL << q; if (r8 == 3 || r8 == 5) m5 |= 1ULL << q; }
    int vf[MAXF];
    for (int j = 1; j < n; j++) {
        u64 l = B[j];
        int v[3] = {0, 0, 0};
        for (int m = 0; m < nf; m++) {
            vf[m] = 0;
            if (!used[m]) continue;
            for (int i = 0; i < ff[m].n; i++) if (ff[m].p[i] == l) { vf[m] = ff[m].e[i]; break; }
        }
        for (int t = 0; t < 3; t++) {
            int s = (l == 2) ? C->c2[t] : 0;
            for (int m = 0; m < nf; m++) s += C->e[t][m] * vf[m];
            v[t] = s;
        }
        int mn = v[0]; if (v[1] < mn) mn = v[1]; if (v[2] < mn) mn = v[2];
        int nm_ = mn & 1, idx = 0, gap = 0;
        for (int i = 0; i < 3; i++) if (v[i] > mn) { idx = i + 1; gap = v[i] - mn; }
        /* unit residues of D1..D3 modulo M (l = 2, 3); for l >= 5 only their Legendre symbols are used,
           so U[t] is computed up to squares: sgn * 2^(c2 mod 2) * prod of the forms with odd exponent */
        u64 M = (l == 2) ? T2_K : (l == 3) ? T3_K : l;
        u64 U[3];
        if (l >= 5) {
            u64 fm[MAXF];
            for (int m = 0; m < nf; m++) {
                if (!((C->e[0][m] | C->e[1][m] | C->e[2][m]) & 1)) continue;
                i64 F = C->F[m]; int neg = F < 0; u64 a = neg ? (u64)(-F) : (u64)F;
                for (int k = 0; k < vf[m]; k++) a = udiv(a, l);
                u64 am = umod(a, l); fm[m] = neg ? l - am : am;
            }
            for (int t = 0; t < 3; t++) {
                u64 r = C->sgn[t] > 0 ? 1 : l - 1;
                if (C->c2[t] & 1) r = mulm(r, 2, l);
                for (int m = 0; m < nf; m++) if (C->e[t][m] & 1) r = mulm(r, fm[m], l);
                U[t] = r;
            }
        } else {
            u64 fm[MAXF];
            for (int m = 0; m < nf; m++) {
                if (!used[m]) continue;
                i64 F = C->F[m]; int neg = F < 0; u64 a = neg ? (u64)(-F) : (u64)F;
                for (int k = 0; k < vf[m]; k++) a = udiv(a, l);
                u64 am = umod(a, M); fm[m] = neg ? (M - am) % M : am;
            }
            for (int t = 0; t < 3; t++) {
                u64 r = C->sgn[t] > 0 ? 1 : M - 1;
                if (l != 2) for (int k = 0; k < C->c2[t]; k++) r = mulm(r, 2, M);
                for (int m = 0; m < nf; m++) for (int k = 0; k < C->e[t][m]; k++) r = mulm(r, fm[m], M);
                U[t] = r;
            }
        }
        int masks[3], nm, hd;
        if (l == 2 || l == 3) {
            int cap = (l == 2) ? T2_CAP : T3_CAP;
            int i1 = (idx == 1) ? 1 : 0;
            u64 a = U[i1];
            u64 tt = idx ? U[idx - 1] : U[(i1 == 0) ? 1 : 2];
            int bc = gap < cap ? gap : cap + (gap - cap) % 2;
            const uint8_t *ent;
            if (l == 2) { ent = T->tab2 + 3 * ((((nm_ * 4 + idx) * (T2_CAP + 2) + bc) * T2_K + a) * T2_K + tt); hd = 3; }
            else        { ent = T->tab3 + 2 * ((((nm_ * 4 + idx) * (T3_CAP + 2) + bc) * T3_K + a) * T3_K + tt); hd = 2; }
            if (ent[0] == 255) return -11;
            nm = hd; for (int k = 0; k < hd; k++) masks[k] = ent[k];
        } else {
            hd = 2;
            if (nm_ != 0 || idx == 0) return -12;
            int i, jj, k; u64 eps;
            if (idx == 1)      { i = 0; jj = 1; k = 2; eps = (l - U[1]) % l; }
            else if (idx == 2) { i = 0; jj = 2; k = 1; eps = (l - U[0]) % l; }
            else               { i = 1; jj = 2; k = 0; eps = U[0]; }
            int split = legendre(eps, l) == 1;
            int tri1[3] = {0, 0, 0}, tri2[3] = {0, 0, 0};
            if (split) { tri1[i] = tri1[jj] = 1; tri2[i] = tri2[jj] = 2; }
            else if (!(gap & 1)) { tri1[0] = 2; tri1[2] = 2; tri2[1] = 2; tri2[2] = 2; }
            else {
                u64 aa = (l - U[idx - 1]) % l;
                tri1[i] = tri1[jj] = 2;
                int la = legendre(aa, l) == -1, le = 1;
                tri2[i] = 1 | ((la ^ le) << 1); tri2[jj] = 1 | (la << 1); tri2[k] = le << 1;
            }
            int g1 = tri1[0] | (tri1[1] << 2), g2 = tri2[0] | (tri2[1] << 2);
            int out[2]; annih4(g1, g2, out); masks[0] = out[0]; masks[1] = out[1]; nm = 2;
        }
        mk[j] = (u32)masks[0] | ((u32)masks[1] << 6) | ((u32)(nm == 3 ? masks[2] : 0) << 12);
        int nl = 0; if (l != 2) for (int kk = 0; kk < nm; kk++) if (masks[kk] & 0xA) nl = 1;
        meta[j] = (uint8_t)(nm | (hd << 2) | (nl << 4));
    }
    /* rows into an echelon basis; primes whose rows need no Legendre symbol first; stop as soon as the
       upper bound 2n - rank on dim Sel2 falls below stop_dim (then the returned value is < stop_dim, not exact) */
#if defined(STAGE) && STAGE == 2
    return n + (int)(mk[1] & 1);
#endif
    row_t piv[2 * MAXB]; u64 pmask[2] = {0, 0}; int rank = 0;
    row_t *P = Pout; int *pv = pvout;           /* optional echelon output (host), NULL on the device */
    u64 Lm[MAXB]; u64 lmdone = 0;                 /* Legendre masks of processed primes, reused by reciprocity */
    u64 m3 = 0;                                   /* bit q: B[q] = 3 mod 4 (odd primes) */
    for (int q = 2; q < n; q++) if ((B[q] & 3) == 3) m3 |= 1ULL << q;
    insert_row_t(piv, pmask, &rank, R0, P, pv);
    for (int pass = 0; pass < 2; pass++) for (int j = 1; j < n; j++) {
        int nl = (meta[j] >> 4) & 1; if (nl != pass) continue;
        int nm = meta[j] & 3, hd = (meta[j] >> 2) & 3;
        u64 l = B[j];
        u64 cm[3];
        if (l == 2) { cm[0] = 2; cm[1] = m1; cm[2] = m5; }
        else {
            cm[0] = 1ULL << j; cm[2] = 0;
            u64 lm = 0;
            if (nl) {
                if ((l & 3) == 3) lm |= 1;
                { int r8 = (int)(l & 7); if (r8 == 3 || r8 == 5) lm |= 2; }
                int l3 = (l & 3) == 3;
                for (int q = 2; q < n; q++) {
                    if (q == j) continue;
                    int bit;
                    if ((lmdone >> q) & 1) bit = (int)((Lm[q] >> j) & 1) ^ (l3 & (int)((m3 >> q) & 1));   /* (b|l) = (l|b) * (-1)^.. */
                    else bit = nonres(B[q], l, T);
                    if (bit) lm |= 1ULL << q;
                }
                Lm[j] = lm; lmdone |= 1ULL << j;
            }
            cm[1] = lm;
        }
        int lowm = (1 << hd) - 1;
        for (int kk = 0; kk < nm; kk++) {
            int mask = (int)((mk[j] >> (6 * kk)) & 63);
            int alo = mask & lowm, ahi = (mask >> hd) & lowm;
            u64 A = 0, Bm = 0;
            for (int bb = 0; bb < hd; bb++) { if ((alo >> bb) & 1) A ^= cm[bb]; if ((ahi >> bb) & 1) Bm ^= cm[bb]; }
            row_t x;
            x.lo = A | (Bm << n);
            x.hi = Bm >> (64 - n);          /* 2 <= n <= MAXB < 64 */
#if defined(STAGE) && STAGE == 3
            rank += (int)((x.lo ^ x.hi) & 1);
#else
            insert_row_t(piv, pmask, &rank, x, P, pv);
#endif
        }
        if (2 * n - rank < stop_dim) return 2 * n - rank;
    }
    if (Bout) { for (int q = 0; q < n; q++) Bout[q] = B[q]; *nout = n; *rankout = rank; }
    return 2 * n - rank;
}
HDN int sel2_dim(const curve_t *C, const fact_t *ff, const seltab_t *T, int stop_dim) {
    return sel2_dim_ex(C, ff, T, stop_dim, (u64 *)0, (int *)0, (row_t *)0, (int *)0, (int *)0);
}


/* ---------------- 3-isogeny Selmer bound, z2z6 (closed form validated against oracle3.gp) ----------------
 * E has the rational 3-torsion point T; phi: E -> E/<T>, phihat dual (kernel mu_3).  Sel^phihat in (Q^x)/(Q^x)^3:
 *   primes of the width-6 forms q, p-5q, p-q (l >= 5): no condition;  primes of the width-2 forms p-3q, p+3q, p-9q
 *   (l >= 5): d must be a local cube (v_l = 0 mod 3; cubic character trivial if l = 1 mod 3);
 *   l = 2: no condition iff q even or p even or p = q or 5q mod 16, else v_2(d) = 0 mod 3;
 *   l = 3: no condition iff 3 does not divide p, else d is a 3-adic cube (v_3 = 0 mod 3, unit = +-1 mod 9).
 * Sel^phi from Greenberg-Wiles: s_phi = s_phihat + corr, corr = -t2 + (1-t3) - #S6 + #{l in S2, l = 1 mod 3}.
 * B3 = s_phi + s_phihat - 1 >= rank. */
/* Barrett modular power for l < 2^32: x mod l = x - mulhi(x, minv) l (+ one correction), minv = floor((2^64-1)/l) */
HD u64 powmod_b32(u64 a, u64 e, u64 l) {
    u64 minv = (~0ULL) / l, r = 1, x = umod(a, l);
    while (e) {
        if (e & 1) { u64 t = r * x; t -= mulhi64(t, minv) * l; if (t >= l) t -= l; r = t; }
        { u64 t = x * x; t -= mulhi64(t, minv) * l; if (t >= l) t -= l; x = t; }
        e >>= 1;
    }
    return r;
}
HD int cubechar(u64 b, u64 l) {              /* l = 1 mod 3 prime, b not divisible by l: 0,1,2 */
    u64 w = (l < (1ULL << 32)) ? powmod_b32(b, (l - 1) / 3, l) : powmod_any(b % l, (l - 1) / 3, l);
    if (w == 1) return 0;
    u64 w2 = l - 1 - w;                       /* the other primitive cube root of unity */
    return (w < w2) ? 1 : 2;
}
HD int cls9(u64 b) { int u = (int)(b % 9); return (u == 2 || u == 7) ? 1 : (u == 4 || u == 5) ? 2 : 0; }
#define MAXB3 48
/* z12 (curve with Z/6 x Z/2 in the class; same B3 for the whole 2-power orbit):
 *   full image: primes l>=5 of p, p-q, u=2p-q, 2p^2-2pq+q^2 (cusp widths 12,12,6,3);
 *   trivial image: primes l>=5 of q, 3p^2-3pq+q^2, 6p^2-6pq+q^2 (widths 2,4,1);
 *   l = 2: no condition unless 4 | q;  l = 3: no condition unless 3 | q (then d must be a 3-adic cube). */
HDN int b3_fam(int fam, i64 p, i64 q, const fact_t *ff) {
    int w6[4], w2[3], n6f, n2f, t2, t3;
    if (fam == FAM_Z2Z6) {      /* family_forms(z2z6): 0 q, 1 p-3q, 2 p-5q, 3 p+3q, 4 p-9q, 5 p-q */
        w6[0] = 0; w6[1] = 2; w6[2] = 5; n6f = 3; w2[0] = 1; w2[1] = 3; w2[2] = 4; n2f = 3;
        t2 = ((q & 1) == 0 || (p & 1) == 0 || ((p - q) & 15) == 0 || ((p - 5 * q) & 15) == 0) ? 1 : 0;
        t3 = (p % 3 == 0) ? 0 : 2;
    } else if (fam == FAM_Z12) { /* family_forms(z12): 0 q, 1 u, 2 A2, 3 A6, 4 p, 5 p-q, 6 A3 */
        w6[0] = 1; w6[1] = 2; w6[2] = 4; w6[3] = 5; n6f = 4; w2[0] = 0; w2[1] = 3; w2[2] = 6; n2f = 3;
        t2 = ((q & 3) == 0) ? 0 : 1;
        t3 = (q % 3 == 0) ? 0 : 2;
    } else return -21;
    u64 B[MAXB3]; int n = 2, n6 = 0; B[0] = 2; B[1] = 3;
    for (int k = 0; k < n6f; k++) for (int i = 0; i < ff[w6[k]].n; i++) {
        u64 l = ff[w6[k]].p[i]; if (l < 5) continue;
        int seen = 0; for (int j = 2; j < n; j++) if (B[j] == l) seen = 1;
        if (!seen) { if (n >= MAXB3) return -20; B[n++] = l; n6++; }
    }
    int s2start = n;
    for (int k = 0; k < n2f; k++) for (int i = 0; i < ff[w2[k]].n; i++) {
        u64 l = ff[w2[k]].p[i]; if (l < 5) continue;
        int seen = 0; for (int j = 2; j < n; j++) if (B[j] == l) seen = 1;
        if (!seen) { if (n >= MAXB3) return -20; B[n++] = l; }
    }
    uint8_t M[2 * MAXB3 + 4][MAXB3]; int nr = 0;
    int corr = -t2 + (1 - t3) - n6;
    #define NEWROW() do { for (int c_ = 0; c_ < n; c_++) M[nr][c_] = 0; } while (0)
    if (!t2) { NEWROW(); M[nr][0] = 1; nr++; }
    if (!t3) { NEWROW(); M[nr][1] = 1; nr++; NEWROW(); for (int c = 0; c < n; c++) if (c != 1) M[nr][c] = (uint8_t)cls9(B[c]); nr++; }
    for (int j = s2start; j < n; j++) {
        u64 l = B[j];
        NEWROW(); M[nr][j] = 1; nr++;
        if (l % 3 == 1) { corr++; NEWROW(); for (int c = 0; c < n; c++) if (c != j) M[nr][c] = (uint8_t)cubechar(B[c], l); nr++; }
    }
    #undef NEWROW
    /* rank mod 3 */
    int rank = 0;
    for (int c = 0; c < n && rank < nr; c++) {
        int piv = -1; for (int r = rank; r < nr; r++) if (M[r][c]) { piv = r; break; }
        if (piv < 0) continue;
        if (piv != rank) for (int k = 0; k < n; k++) { uint8_t t = M[piv][k]; M[piv][k] = M[rank][k]; M[rank][k] = t; }
        int inv = M[rank][c];                 /* 1 or 2; inverse mod 3 is itself */
        for (int k = 0; k < n; k++) M[rank][k] = (uint8_t)((M[rank][k] * inv) % 3);
        for (int r = 0; r < nr; r++) if (r != rank && M[r][c]) {
            int f = M[r][c];
            for (int k = 0; k < n; k++) M[r][k] = (uint8_t)((M[r][k] + 3 * 3 - f * M[rank][k]) % 3);
        }
        rank++;
    }
    int sh = n - rank;
    return 2 * sh + corr - 1;
}


/* B3 with early exit (target > 0): exact when the result is >= target, otherwise some value < target.
 * After the valuation rows (one per width-2 prime, e_2 if t2 = 0, e_3 if t3 = 0) only the columns
 * R = {2 if t2} u {3 if t3 != 0} u {width-6 primes} matter; UB = |R|-based bound = n6 + n21 + t2, and every
 * independent extra row (3-adic unit row, cubic-character rows of width-2 primes l = 1 mod 3) lowers B3 by 2. */
HDN int b3_fast(int fam, i64 p, i64 q, const fact_t *ff, int target, int par) {
    int w6[4], w2[3], n6f, n2f, t2, t3;
    if (fam == FAM_Z2Z6) {
        w6[0] = 0; w6[1] = 2; w6[2] = 5; n6f = 3; w2[0] = 1; w2[1] = 3; w2[2] = 4; n2f = 3;
        t2 = ((q & 1) == 0 || (p & 1) == 0 || ((p - q) & 15) == 0 || ((p - 5 * q) & 15) == 0) ? 1 : 0;
        t3 = (p % 3 == 0) ? 0 : 2;
    } else if (fam == FAM_Z12) {
        w6[0] = 1; w6[1] = 2; w6[2] = 4; w6[3] = 5; n6f = 4; w2[0] = 0; w2[1] = 3; w2[2] = 6; n2f = 3;
        t2 = ((q & 3) == 0) ? 0 : 1;
        t3 = (q % 3 == 0) ? 0 : 2;
    } else return -21;
    u64 Rc[MAXB3]; int nR = 0;                 /* the surviving columns */
    if (t2) Rc[nR++] = 2;
    if (t3) Rc[nR++] = 3;
    int n6 = 0;
    for (int k = 0; k < n6f; k++) for (int i = 0; i < ff[w6[k]].n; i++) {
        u64 l = ff[w6[k]].p[i]; if (l < 5) continue;
        int seen = 0; for (int j = 0; j < nR; j++) if (Rc[j] == l) seen = 1;
        if (!seen) { if (nR >= MAXB3) return -20; Rc[nR++] = l; n6++; }
    }
    u64 L1[MAXB3]; int n21 = 0;                /* width-2 primes = 1 mod 3 (character rows) */
    for (int k = 0; k < n2f; k++) for (int i = 0; i < ff[w2[k]].n; i++) {
        u64 l = ff[w2[k]].p[i]; if (l < 5 || l % 3 != 1) continue;
        int seen = 0; for (int j = 0; j < n21; j++) if (L1[j] == l) seen = 1;
        if (!seen) { if (n21 >= MAXB3) return -20; L1[n21++] = l; }
    }
    int corr = -t2 + (1 - t3) - n6 + n21;
    int rank = 0;
    #define B3CUR (2 * (nR - rank) + corr - 1)
#if defined(STAGE) && STAGE == 10
    return -1 - (B3CUR & 0);
#endif
    if (target > 0 && B3CUR < target) return B3CUR;
    /* B3 = UB mod 2 (each independent row lowers it by 2), and B2 = B3 = rank parity mod 2: wrong parity -> reject */
    if (target > 0 && par >= 0 && ((B3CUR ^ par) & 1)) return target - 1;
    uint8_t P[MAXB3][MAXB3]; int pc[MAXB3];
    uint8_t x[MAXB3];
    for (int r = -1; r < n21; r++) {
        if (r < 0) { if (t3) continue; for (int c = 0; c < nR; c++) x[c] = (uint8_t)cls9(Rc[c]); }   /* 3-adic unit row */
        else if (rank == 0 && target > 0) {
            /* first row: its first nonzero entry already proves rank 1; exit before computing the rest if that suffices */
            int c = 0;
            for (; c < nR; c++) { x[c] = (uint8_t)cubechar(Rc[c], L1[r]); if (x[c]) break; }
            if (c < nR && 2 * (nR - 1) + corr - 1 < target) return 2 * (nR - 1) + corr - 1;
            for (c = c + 1; c < nR; c++) x[c] = (uint8_t)cubechar(Rc[c], L1[r]);
        }
        else for (int c = 0; c < nR; c++) x[c] = (uint8_t)cubechar(Rc[c], L1[r]);
        for (int k = 0; k < rank; k++) { int f = x[pc[k]]; if (f) for (int c = 0; c < nR; c++) x[c] = (uint8_t)((x[c] + 3 * 3 - f * P[k][c]) % 3); }
        int c0 = -1; for (int c = 0; c < nR; c++) if (x[c]) { c0 = c; break; }
        if (c0 < 0) continue;
        int inv = x[c0];                        /* 1 or 2, self-inverse mod 3 */
        for (int c = 0; c < nR; c++) P[rank][c] = (uint8_t)((x[c] * inv) % 3);
        pc[rank] = c0; rank++;
        if (target > 0 && B3CUR < target) return B3CUR;
    }
    return B3CUR;
    #undef B3CUR
}




/* z8: predicted rank parity from the factorisation (root number w = -w_2 prod_{l odd} w_l):
 * odd primes of p, p-q, 2p-q: split (w_l = -1); of q: split iff l = 1 mod 4; of Q = 8p^2-8pq+q^2: split iff l = 1 mod 8;
 * w_2 from tabw2 by the tab2 key of E1 (validated: genw2.gp). Returns 0/1 (rank parity), or -1 if unknown. */
HDN int z8_parity(const curve_t *C, const fact_t *ff, const uint8_t *tabw2) {
    int neg = 1;                                  /* w_inf = -1 */
    for (int m = 0; m < 5; m++) for (int i = 0; i < ff[m].n; i++) {
        u64 l = ff[m].p[i]; if (l == 2) continue;
        int split = (m <= 2) ? 1 : (m == 3) ? ((l & 3) == 1) : ((l & 7) == 1);
        if (split) neg ^= 1;
    }
    /* 2-adic key of E1 (same as sel2_dim's l = 2 branch) */
    int vf[MAXF], v[3];
    for (int m = 0; m < C->nf; m++) { vf[m] = 0; for (int i = 0; i < ff[m].n; i++) if (ff[m].p[i] == 2) vf[m] = ff[m].e[i]; }
    for (int t = 0; t < 3; t++) { int s_ = C->c2[t]; for (int m = 0; m < C->nf; m++) s_ += C->e[t][m] * vf[m]; v[t] = s_; }
    int mn = v[0]; if (v[1] < mn) mn = v[1]; if (v[2] < mn) mn = v[2];
    int nm_ = mn & 1, idx = 0, gap = 0;
    for (int i = 0; i < 3; i++) if (v[i] > mn) { idx = i + 1; gap = v[i] - mn; }
    if (idx == 0) return -1;
    u64 U[3];
    for (int t = 0; t < 3; t++) {
        u64 r = C->sgn[t] > 0 ? 1 : T2_K - 1;
        for (int m = 0; m < C->nf; m++) { int e = C->e[t][m]; if (!e) continue;
            i64 F = C->F[m]; int ng = F < 0; u64 a = ng ? (u64)(-F) : (u64)F; for (int k = 0; k < vf[m]; k++) a >>= 1;
            u64 am = a % T2_K; if (ng) am = (T2_K - am) % T2_K; for (int k = 0; k < e; k++) r = (r * am) % T2_K; }
        U[t] = r;
    }
    int i1 = (idx == 1) ? 1 : 0; u64 a = U[i1], tt = U[idx - 1];
    int bc = gap < T2_CAP ? gap : T2_CAP + (gap - T2_CAP) % 2;
    int w2 = tabw2[(((nm_ * 4 + idx) * (T2_CAP + 2) + bc) * T2_K + a) * T2_K + tt];
    if (w2 == 0) return -1;
    if (w2 == 2) neg ^= 1;
    return neg;                                   /* w = -1 <=> odd rank parity */
}

/* phase-1 prefilter for the two-phase GPU kernel (B3-first families z2z6, z12, with target > 0):
 * returns 1 if the pair needs the full exact computation, 0 if it is rejected (bound < target or wrong parity),
 * SEL_SINGULAR / SEL_FACTFAIL / SEL_ERROR codes otherwise. Uses only the factorisation (no symbols). */
HDN int prefilter_b3ub(int fam, i64 p, i64 q, int all, const u32 *spf, u64 lim, int target, const uint8_t *tabw2) {
    i64 F[MAXF];
    {   i64 ap = p < 0 ? -p : p, aq = q < 0 ? -q : q, h = ap > aq ? ap : aq;
        if (h > (fam == FAM_Z2Z6 ? (i64)900000000000000000LL : (i64)700000000LL)) return SEL_ERROR; }
    int nf = family_forms(fam, p, q, F);
    if (nf <= 0) return SEL_SINGULAR;
    fact_t ff[MAXF];
    for (int m = 0; m < nf; m++) { u64 a = F[m] < 0 ? (u64)(-F[m]) : (u64)F[m]; if (factor64(a, &ff[m], spf, lim)) return SEL_FACTFAIL; }
    if (fam == FAM_Z8) {                          /* z8: parity cut only */
        if (!(all & 4) || !tabw2) return 1;
        curve_t cv[3]; family_curves(fam, F, nf, cv);
        int pr = z8_parity(&cv[0], ff, tabw2);
        return (pr >= 0 && pr != ((all >> 3) & 1)) ? 0 : 1;
    }
    int w6[4], w2[3], n6f, n2f, t2;
    if (fam == FAM_Z2Z6) { w6[0] = 0; w6[1] = 2; w6[2] = 5; n6f = 3; w2[0] = 1; w2[1] = 3; w2[2] = 4; n2f = 3;
        t2 = ((q & 1) == 0 || (p & 1) == 0 || ((p - q) & 15) == 0 || ((p - 5 * q) & 15) == 0) ? 1 : 0; }
    else { w6[0] = 1; w6[1] = 2; w6[2] = 4; w6[3] = 5; n6f = 4; w2[0] = 0; w2[1] = 3; w2[2] = 6; n2f = 3; t2 = ((q & 3) == 0) ? 0 : 1; }
    /* UB = n6 + n21 + t2 (distinct primes >= 5; primes >= 5 of different forms of one class never coincide here,
       but dedupe anyway) */
    u64 seen[MAXB3]; int ns = 0, n6 = 0, n21 = 0;
    for (int k = 0; k < n6f; k++) for (int i = 0; i < ff[w6[k]].n; i++) { u64 l = ff[w6[k]].p[i]; if (l < 5) continue;
        int d = 0; for (int j = 0; j < ns; j++) if (seen[j] == l) d = 1; if (!d && ns < MAXB3) { seen[ns++] = l; n6++; } }
    ns = 0;
    for (int k = 0; k < n2f; k++) for (int i = 0; i < ff[w2[k]].n; i++) { u64 l = ff[w2[k]].p[i]; if (l < 5 || l % 3 != 1) continue;
        int d = 0; for (int j = 0; j < ns; j++) if (seen[j] == l) d = 1; if (!d && ns < MAXB3) { seen[ns++] = l; n21++; } }
    int ub = n6 + n21 + t2;
    if (ub < target) return 0;
    if ((all & 4) && ((ub ^ (all >> 3)) & 1)) return 0;
    return 1;
}


/* ---------------- z10: 2-isogeny descent pair, closed form (validated against oracle.gp selphi) ----------------
 * E: y^2 = x(x^2 + A x + B), A = -(2p^2-2pq+q^2)(4p^4-12p^3q+6p^2q^2+2pq^3-q^4), B = 16 p^5 (p-q)^5 P2,
 * E': y^2 = x(x^2 - 2A x + B'), B' = A^2 - 4B = q^5 u^5 F4.  W1_v, W2_v = images of E(Q_v), E'(Q_v) under x -> x mod squares.
 *   odd l != 5, l | B (primes of p, p-q, P2): split ((A|l) = 1) or v_l(B) odd: W1 = all, W2 = 0; else W1 = W2 = units.
 *   odd l != 5, l | B' (q, u, F4): split ((-2A|l) = 1) or v_l(B') odd: W1 = 0, W2 = all; else W1 = W2 = units.
 *   l = 2: 2 | q: W1 = 0, W2 = all; otherwise W1 = all, W2 = 0.
 *   l = 5: 5 | p(p-q): W1 = all, W2 = 0; 5 | q u: W1 = 0, W2 = all; 5 | 2p^2-2pq+q^2: W1 = W2 = units; 5 | P2: W1 = W2 = <5*2>.
 *   infinity: W = all iff the cubic has a negative root (B < 0, or B' > 0 and A > 0 for E; same with (-2A, B') for E').
 * returns s1 + s2 - 2 (rank bound), s1/s2 via *sd; negative on error. */
HD u64 z10_A_mod(i64 p, i64 q, u64 l) {
    u64 pm = (u64)((p % (i64)l + (i64)l) % (i64)l), qm = (u64)((q % (i64)l + (i64)l) % (i64)l);
    #define MM(a, b) mulmod_any((a), (b), l)
    u64 p2 = MM(pm, pm), q2 = MM(qm, qm), pq = MM(pm, qm);
    u64 f1 = (MM(2, p2) + l - MM(2, pq) + q2) % l;                       /* 2p^2 - 2pq + q^2 */
    u64 p3 = MM(p2, pm), p4 = MM(p3, pm), q3 = MM(q2, qm), q4 = MM(q3, qm);
    u64 k4 = (MM(4, p4) + l - MM(12, MM(p3, qm)) + MM(6, MM(p2, q2)) + MM(2, MM(pm, q3)) + l - q4) % l;
    #undef MM
    return (l - mulmod_any(f1, k4, l)) % l;                              /* A = -f1 k4 */
}
HD u64 z10b_A_mod(i64 p, i64 q, u64 l) {
    u64 pm = (u64)((p % (i64)l + (i64)l) % (i64)l), qm = (u64)((q % (i64)l + (i64)l) % (i64)l);
    #define MM(a, b) mulmod_any((a), (b), l)
    u64 p2 = MM(pm, pm), q2 = MM(qm, qm), pq = MM(pm, qm);
    u64 f1 = (MM(2, p2) + l - MM(2, pq) + q2) % l;
    u64 p3 = MM(p2, pm), p4 = MM(p3, pm), q3 = MM(q2, qm), q4 = MM(q3, qm);
    u64 k4 = (MM(4, p4) + MM(36 % l, MM(p3, qm)) + l - MM(66 % l, MM(p2, q2)) + MM(26 % l, MM(pm, q3)) + l - q4) % l;
    #undef MM
    return (l - mulmod_any(f1, k4, l)) % l;                              /* A5 = -f1 k4b */
}
/* v5 = 1: the 5-isogenous pair (E5, E5') = (E/<2P>, E'/<2P'>): A5 = -(2p^2-2pq+q^2)(4p^4+36p^3q-66p^2q^2+26pq^3-q^4),
 * B5 = 16 p (p-q) P2^5, B5' = q u F4^5; the local rules are the same as for (E, E') (brute-force oracle, work/selmer/rules10b.out). */
HDN int bphi_z10_gen(int v5, i64 p, i64 q, const fact_t *ff, const seltab_t *T, int *sd, u64 *Bout, int *nout, u64 *R1out, int *n1out, u64 *R2out, int *n2out) {
    /* forms: 0 p, 1 p-q, 2 P2, 3 q, 4 u, 5 F4 ; class of a prime: 1 = divides B, 2 = divides B' */
    u64 B_[MAXB]; int cl[MAXB], vv[MAXB], n = 2;
    B_[0] = (u64)-1; B_[1] = 2; cl[0] = cl[1] = 0;
    static const int exE[6] = {5, 5, 1, 5, 5, 1}, exE5[6] = {1, 1, 5, 1, 1, 5}; const int *ex = v5 ? exE5 : exE;
    for (int m = 0; m < 6; m++) for (int i = 0; i < ff[m].n; i++) {
        u64 l = ff[m].p[i]; if (l == 2) continue;
        int j = -1; for (int k = 2; k < n; k++) if (B_[k] == l) j = k;
        if (j < 0) { if (n >= MAXB) return -10; j = n++; B_[j] = l; cl[j] = 0; vv[j] = 0; }
        cl[j] |= (m < 3) ? 1 : 2; vv[j] += ex[m] * ff[m].e[i];
    }
    /* exact A, B, B' signs (64-bit overflow-free via signs of the forms) */
    i64 f1 = 2 * p * p - 2 * p * q + q * q;                               /* > 0 always */
    /* sign(K4) needs care: evaluate in floating point is unsafe; use 128-bit integers */
    __int128 P = p, Q = q;
    __int128 K4 = v5 ? 4 * P * P * P * P + 36 * P * P * P * Q - 66 * P * P * Q * Q + 26 * P * Q * Q * Q - Q * Q * Q * Q
                     : 4 * P * P * P * P - 12 * P * P * P * Q + 6 * P * P * Q * Q + 2 * P * Q * Q * Q - Q * Q * Q * Q;
    int sA = -(f1 > 0 ? 1 : -1) * (K4 > 0 ? 1 : (K4 < 0 ? -1 : 0));
    i64 P2 = p * p - 3 * p * q + q * q, F4 = 4 * p * p - 2 * p * q - q * q;
    int sB = ((p < 0) ^ (p - q < 0) ^ (P2 < 0)) ? -1 : 1;                /* B = 16 p^5 (p-q)^5 P2 */
    int sB2 = ((q < 0) ^ (2 * p - q < 0) ^ (F4 < 0)) ? -1 : 1;           /* B' = q^5 u^5 F4 */
    if (sA == 0) return -13;
    int neg1 = (sB < 0) || (sB2 > 0 && sA > 0);                           /* E has a real point with x < 0 */
    int neg2 = (sB2 < 0) || (sB > 0 && sA < 0);                           /* E': (A',B'') = (-2A, A'^2-4B' = 16 B) */
    /* rows over n columns (bit q = basis element q) */
    u64 R1[3 * MAXB + 8], R2[3 * MAXB + 8]; int n1 = 0, n2 = 0;
    u64 m1 = 1, m5 = 0;                                                    /* 2-adic class masks of the basis */
    for (int k = 2; k < n; k++) { int r8 = (int)(B_[k] & 7); if ((r8 & 3) == 3) m1 |= 1ULL << k; if (r8 == 3 || r8 == 5) m5 |= 1ULL << k; }
    if (!neg1) R1[n1++] = 1;                                              /* d > 0 */
    if (!neg2) R2[n2++] = 1;
    if ((q & 1) == 0) { R1[n1++] = 2; R1[n1++] = m1; R1[n1++] = m5; }     /* W1 = 0 at 2 */
    else              { R2[n2++] = 2; R2[n2++] = m1; R2[n2++] = m5; }
    for (int j = 2; j < n; j++) {
        u64 l = B_[j];
        u64 lm = 0;                                                        /* Legendre mask: bit k iff B_[k] non-residue mod l */
        if ((l & 3) == 3) lm |= 1;
        { int r8 = (int)(l & 7); if (r8 == 3 || r8 == 5) lm |= 2; }
        for (int k = 2; k < n; k++) if (k != j && nonres(B_[k], l, T)) lm |= 1ULL << k;
        u64 vrow = 1ULL << j;
        if (l == 5) {
            u64 d5 = (u64)(((p % 5) + 5) % 5), q5 = (u64)(((q % 5) + 5) % 5);
            int c;                                                         /* which form 5 divides */
            if (d5 == 0 || (d5 + 5 - q5) % 5 == 0) c = 1;                  /* p, p-q */
            else if (q5 == 0 || (2 * d5 + 5 - q5) % 5 == 0) c = 2;         /* q, u */
            else { u64 p2 = (d5 * d5 + 5 * 5 * 3 - 3 * d5 * q5 + q5 * q5) % 5; c = (p2 == 0) ? 4 : 3; }   /* P2 (p = -q mod 5; A7 also vanishes) or the A-forms only (p = 2q) */
            if (c == 1) { R2[n2++] = vrow; R2[n2++] = lm; }
            else if (c == 2) { R1[n1++] = vrow; R1[n1++] = lm; }
            else if (c == 3) { R1[n1++] = vrow; R2[n2++] = vrow; }
            else { R1[n1++] = vrow ^ lm; R2[n2++] = vrow ^ lm; }          /* annihilator of <(v=1, leg=1)> */
            continue;
        }
        if (cl[j] == 3) return -14;                                       /* l | B and l | B': not in this family for l != 5 */
        u64 Am = v5 ? z10b_A_mod(p, q, l) : z10_A_mod(p, q, l);
        if (cl[j] == 1) {
            int split = legendre(Am, l) == 1;
            if (split || (vv[j] & 1)) { R2[n2++] = vrow; R2[n2++] = lm; }
            else { R1[n1++] = vrow; R2[n2++] = vrow; }
        } else {
            int split = legendre(mulmod_any(l - Am, 2, l), l) == 1;       /* (-2A | l) */
            if (split || (vv[j] & 1)) { R1[n1++] = vrow; R1[n1++] = lm; }
            else { R1[n1++] = vrow; R2[n2++] = vrow; }
        }
    }
    /* ranks over F2 (n <= 56 columns: one word) */
    int rk[2] = {0, 0};
    for (int g = 0; g < 2; g++) {
        u64 *R = g ? R2 : R1; int nr = g ? n2 : n1; u64 piv[64]; u64 pm = 0;
        for (int i = 0; i < nr; i++) { u64 y = R[i];
            while (y) { int c = ctz64(y); if (!((pm >> c) & 1)) { piv[c] = y; pm |= 1ULL << c; rk[g]++; break; } y ^= piv[c]; } }
    }
    int s1 = n - rk[0], s2 = n - rk[1];
    if (sd) { sd[0] = s1; sd[1] = s2; }
    if (Bout) { for (int k = 0; k < n; k++) Bout[k] = B_[k]; *nout = n;
                for (int i = 0; i < n1; i++) R1out[i] = R1[i]; *n1out = n1; for (int i = 0; i < n2; i++) R2out[i] = R2[i]; *n2out = n2; }
    return s1 + s2 - 2;
}
HDN int bphi_z10_ex(i64 p, i64 q, const fact_t *ff, const seltab_t *T, int *sd, u64 *Bout, int *nout, u64 *R1out, int *n1out, u64 *R2out, int *n2out) {
    return bphi_z10_gen(0, p, q, ff, T, sd, Bout, nout, R1out, n1out, R2out, n2out);
}
HDN int bphi_z10b_ex(i64 p, i64 q, const fact_t *ff, const seltab_t *T, int *sd, u64 *Bout, int *nout, u64 *R1out, int *n1out, u64 *R2out, int *n2out) {
    return bphi_z10_gen(1, p, q, ff, T, sd, Bout, nout, R1out, n1out, R2out, n2out);
}
HDN int bphi_z10(i64 p, i64 q, const fact_t *ff, const seltab_t *T, int *sd) {
    return bphi_z10_ex(p, q, ff, T, sd, (u64 *)0, (int *)0, (u64 *)0, (int *)0, (u64 *)0, (int *)0);
}



/* ---------------- z10: BSD budget estimate bhat = Omega * prod(c_p) / 100 (owner note 10, job F) ----------------
 * Omega = real period omega1 of the family model (search/families/fam.gp, minimal away from 2, 3, 5), by the AGM on the
 * real roots of 4x^3 + b2 x^2 + 2 b4 x + b6 (as PARI's E.omega[1]; ellbsd does not count real components).
 * prod(c_p) from the factorisation: l | p, p-q: 10 v; l | 2p-q, q: 5 v; l | P2 = p^2-3pq+q^2: 2v if split else 2;
 * l | F4 = 4p^2-2pq-q^2: v if split else (v odd ? 1 : 2); split <=> -c6 is a square mod l.  l = 2, 3, 5 and primes dividing two
 * forms: 1 (lower bound).  Returns log10(bhat), or +99 if the period computation fails (so a filter keeps the pair).
 * Validated against PARI ellbsd on 2e4 random pairs (h 1e3..3e5): log10(bhat/b) in [-0.602, 0], 0 overestimates. */
HD double z10_omega1(i64 p, i64 q) {
    /* f(x) = 4x^3 + b2 x^2 + 2 b4 x + b6 = (x - X)(4x^2 + B x + C) with the rational 2-torsion root X = -p^3 (p-q) P2,
     * B = -q^2 u K3 (K3 = 4p^3-4p^2q-2pq^2+q^3), C = q^4 p^3 (p-q) u^2 P2, disc(4x^2+Bx+C) = D = q^5 u^5 F4,
     * f'(X) = G = 4 p^5 (p-q)^5 P2: every root difference from products of forms (no cancellation). */
    double P = (double)p, Q = (double)q, PQ = (double)(p - q), U = (double)(2 * p - q);
    double P2 = P * P - 3 * P * Q + Q * Q, F4 = 4 * P * P - 2 * P * Q - Q * Q, K3 = 4 * P * P * P - 4 * P * P * Q - 2 * P * Q * Q + Q * Q * Q;
    double X = -P * P * P * PQ * P2;
    double B = -Q * Q * U * K3;
    double C = Q * Q * Q * Q * P * P * P * PQ * U * U * P2;
    double D = pow(Q * U, 5) * F4;
    double G = 4 * pow(P * PQ, 5) * P2;
    const double PI = 3.14159265358979323846;
    if (D > 0) {                                          /* three real roots X, r1, r2 */
        double sd = sqrt(D), r1 = (-B - (B >= 0 ? sd : -sd)) / 8, r2 = (r1 != 0) ? C / (4 * r1) : (-B + sd) / 8;
        double hi = r1 > r2 ? r1 : r2, lo = r1 > r2 ? r2 : r1, dhl = sd / 4;          /* hi - lo */
        double sa = X - hi, sb = X - lo;                  /* (X-hi)(X-lo) = G/4: recompute the smaller one */
        if (fabs(sa) < fabs(sb)) sa = (G / 4) / sb; else sb = (G / 4) / sa;
        double d13, d12;                                   /* e1 - e3, e1 - e2 with e1 > e2 > e3 */
        if (sa > 0) { d13 = sb; d12 = sa; }               /* X > hi > lo */
        else if (sb < 0) { d13 = -sa; d12 = dhl; }        /* hi > lo > X: e1 = hi, e3 = X; e1-e2 = hi-lo */
        else { d13 = dhl; d12 = -sa; }                     /* hi > X > lo */
        if (!(d13 > 0 && d12 > 0)) return -1;
        double x = sqrt(d13), y = sqrt(d12);
        for (int it = 0; it < 80 && fabs(x - y) > 1e-15 * x; it++) { double m = (x + y) / 2; y = sqrt(x * y); x = m; }
        return PI / x;
    } else {                                              /* one real root X, complex pair z, conj(z) */
        /* w = X - z: Re w = (8X+B)/8 = -(2p^2-2pq+q^2) K4 / 8 (K4 = 4p^4-12p^3q+6p^2q^2+2pq^3-q^4), Im w = sqrt(-D)/8.
         * AGM(2 sqrt|w|, sqrt(2(|w| + Re w))) with |w| + Re w = Im^2 / (|w| - Re w) when Re w < 0 (no cancellation). */
        double f1 = 2 * P * P - 2 * P * Q + Q * Q, K4 = 4 * P * P * P * P - 12 * P * P * P * Q + 6 * P * P * Q * Q + 2 * P * Q * Q * Q - Q * Q * Q * Q;
        double re = -f1 * K4 / 8, im = sqrt(-D) / 8, aw = hypot(re, im);
        double t = (re >= 0) ? aw + re : im * im / (aw - re);
        double x = 2 * sqrt(aw), y = sqrt(2 * t);
        if (!(y > 0 && x > 0)) return -1;
        for (int it = 0; it < 80 && fabs(x - y) > 1e-15 * x; it++) { double m = (x + y) / 2; y = sqrt(x * y); x = m; }
        return 2 * PI / x;
    }
}
HD u64 z10_c6_mod(i64 p, i64 q, u64 l) {          /* c6 of the family model mod l */
    #define MD(x) ((u64)(((x) % (i64)l + (i64)l) % (i64)l))
    #define M2(a, b) mulmod_any((a), (b), l)
    u64 P = MD(p), Q = MD(q), P2 = M2(P, P), P3 = M2(P2, P), P4 = M2(P3, P), P5 = M2(P4, P), P6 = M2(P5, P), P7 = M2(P6, P);
    u64 Q2 = M2(Q, Q), Q3 = M2(Q2, Q), Q4 = M2(Q3, Q), Q5 = M2(Q4, Q), Q6 = M2(Q5, Q);
    u64 a1 = (M2(2, P3) + l - M2(2, M2(Q, P2)) + l - M2(2, M2(Q2, P)) + Q3) % l;
    u64 a2 = (l - M2(2, M2(Q, P5)) + M2(3, M2(Q2, P4)) + l - M2(Q3, P3)) % l;
    u64 a3 = (l - M2(2, M2(Q2, P7)) + M2(9 % l, M2(Q3, P6)) + l - M2(12 % l, M2(Q4, P5)) + M2(6 % l, M2(Q5, P4)) + l - M2(Q6, P3)) % l;
    u64 b2 = (M2(a1, a1) + M2(4 % l, a2)) % l, b4 = M2(a1, a3), b6 = M2(a3, a3);
    u64 t1 = M2(M2(b2, b2), b2), t2 = M2(36 % l, M2(b2, b4)), t3 = M2(216 % l, b6);
    #undef MD
    #undef M2
    return (l - t1 + t2 + l - t3) % l;
}
HDN double z10_log10_bhat(i64 p, i64 q, const fact_t *ff) {
    double om = z10_omega1(p, q);
    if (!(om > 0)) return 99;                     /* estimator failure: keep the pair (filter-safe) */
    double lc = 0;                                  /* log10 prod c_p */
    if ((q & 1) == 0) {                             /* 2 | q and 2 | 2p-q: c_2 = 5 (v2(q) + v2(2p-q) - 2), the model is not minimal at 2 */
        i64 u = 2 * p - q; int a = 0, b = 0; i64 t = q; while (!(t & 1)) { t >>= 1; a++; } t = u; while (t && !(t & 1)) { t >>= 1; b++; }
        if (a + b - 2 > 0) lc += log10(5.0 * (a + b - 2));
        om *= 2;
    }
    for (int m = 0; m < 6; m++) for (int i = 0; i < ff[m].n; i++) {
        u64 l = ff[m].p[i]; int v = ff[m].e[i];
        if (l == 2 && (q & 1) == 0) continue;
        int other = 0; for (int m2 = 0; m2 < 6; m2++) if (m2 != m) for (int j = 0; j < ff[m2].n; j++) if (ff[m2].p[j] == l) other = 1;
        if (other) continue;
        double c;
        if (m == 0 || m == 1) c = 10.0 * v;
        else if (m == 3 || m == 4) c = 5.0 * v;
        else {
            u64 c6 = (l == 2) ? 0 : z10_c6_mod(p, q, l); int split = c6 && legendre((l - c6) % l, l) == 1;   /* l = 2: treated as non-split (lower bound) */
            if (m == 2) c = split ? 2.0 * v : 2.0;
            else c = split ? (double)v : ((v & 1) ? 1.0 : 2.0);
        }
        lc += log10(c);
    }
    return log10(om) + lc - 2.0;
}

/* log10 bhat from (p,q) alone (factors the six z10 forms); +99 when the factorisation fails (keep) */
HDN double z10_log10_bhat_pq(i64 p, i64 q, const u32 *spf, u64 lim) {
    i64 F[MAXF]; int nf = family_forms(FAM_Z10, p, q, F); if (nf <= 0) return 99;
    fact_t ff[MAXF];
    for (int m = 0; m < nf; m++) if (factor64(F[m] < 0 ? (u64)(-F[m]) : (u64)F[m], &ff[m], spf, lim)) return 99;
    return z10_log10_bhat(p, q, ff);
}

/* ---------------- z9: the two 3-isogeny bounds (chain E -> E/<3P> -> E/<P>), closed form (validated vs oraclen.gp) ----------------
 * which = 0: kernel <3P> on E;  which = 1: kernel <image of P> on E/<3P>.
 * full image (no condition): primes of p, q, p-q (always including 2), and for which = 1 also the primes of A9;
 * trivial image (d a local cube): primes >= 5 of C3, and for which = 0 also of A9;  l = 3 with p = -q mod 3: v_3(d) = 0 mod 3.
 * B = 2 (|R| - rank) + corr - 1, R = columns of the full primes, corr = -#full + #(trivial primes = 1 mod 3). Early exit at target. */
HDN int b3_z9(i64 p, i64 q, const fact_t *ff, int which, int target) {
    u64 Rc[MAXB3]; int nR = 0; u64 L1[MAXB3]; int n21 = 0;
    int nfullf = which ? 4 : 3;
    for (int m = 0; m < 5; m++) for (int i = 0; i < ff[m].n; i++) {
        u64 l = ff[m].p[i];
        if (m == 3 && l == 3) continue;            /* 3 | A9 (p = -q mod 3): units case at 3, not a full column */
        if (m < nfullf) { int seen = 0; for (int j = 0; j < nR; j++) if (Rc[j] == l) seen = 1; if (!seen) { if (nR >= MAXB3) return -20; Rc[nR++] = l; } }
        else if (l >= 5 && l % 3 == 1) { int seen = 0; for (int j = 0; j < n21; j++) if (L1[j] == l) seen = 1; if (!seen) { if (n21 >= MAXB3) return -20; L1[n21++] = l; } }
    }
    /* 2 always divides one of p, q, p-q; 3 is full iff it divides one of them (else 3 | A9, C3: valuation row only) */
    int corr = -nR + n21, rank = 0;
    #define B9CUR (2 * (nR - rank) + corr - 1)
    if (target > 0 && B9CUR < target) return B9CUR;
    uint8_t P[MAXB3][MAXB3]; int pc[MAXB3]; uint8_t x[MAXB3];
    for (int r = 0; r < n21; r++) {
        if (rank == 0 && target > 0) {
            int c = 0;
            for (; c < nR; c++) { x[c] = (uint8_t)cubechar(Rc[c], L1[r]); if (x[c]) break; }
            if (c < nR && 2 * (nR - 1) + corr - 1 < target) return 2 * (nR - 1) + corr - 1;
            for (c = c + 1; c < nR; c++) x[c] = (uint8_t)cubechar(Rc[c], L1[r]);
        } else for (int c = 0; c < nR; c++) x[c] = (uint8_t)cubechar(Rc[c], L1[r]);
        for (int k = 0; k < rank; k++) { int f = x[pc[k]]; if (f) for (int c = 0; c < nR; c++) x[c] = (uint8_t)((x[c] + 9 - f * P[k][c]) % 3); }
        int c0 = -1; for (int c = 0; c < nR; c++) if (x[c]) { c0 = c; break; }
        if (c0 < 0) continue;
        int inv = x[c0]; for (int c = 0; c < nR; c++) P[rank][c] = (uint8_t)((x[c] * inv) % 3);
        pc[rank] = c0; rank++;
        if (target > 0 && B9CUR < target) return B9CUR;
    }
    return B9CUR;
    #undef B9CUR
}


/* ---------------- z8 (E0) and z12 (E): 2-isogeny pair on the curve with ONE rational 2-torsion point ----------------
 * z8:  E0: y^2 = x(x^2 + A x + B), A = K8 = u^4 - 8 w^2, B = 16 w^4 (w = p(p-q)), B' = A^2 - 4B = q^2 u^4 Q.
 * z12: E:  A = A12 (octic), B = 16 p^6 (p-q)^6 A3^2, B' = q^2 u^6 A2^3 A6.   (models checked against the family curves)
 * odd l (incl. 3 for z8): as z10 (l | B: split <=> (A|l) = 1; l | B': split <=> (-2A|l) = 1; split or odd valuation -> full
 *   on one side, trivial on the other; else both units).  l = 2: 2 | p(p-q) -> W(E) all, W(E') 0;  2 | q:
 *   z8:  v2(q)=1,v2(u)=2: <5>, units;  v2(q)=1,v2(u)>=3: 0, all;  v2(q)=2: <3>, <5,-2>;  v2(q)>=3: <-1>, <5,2>
 *   z12: v2(q)=1: 0, all;  v2(q)>=2: <5>, units.
 *   z12, l = 3: 3 | p(p-q): all, 0;  3 | u: 0, all;  3 | q: <3*2>, <3>.
 * (derived from oracle.gp phipair data; validated against it.)  Same outputs as bphi_z10_ex. */
HDN int bphi_one_ex(int fam, i64 p, i64 q, const fact_t *ff, const seltab_t *T, int *sd, u64 *Bout, int *nout, u64 *R1out, int *n1out, u64 *R2out, int *n2out) {
    /* form descriptors: side (1 = divides B, 2 = divides B') and exponent */
    int side[MAXF], ex[MAXF], nf;
    if (fam == FAM_Z8) {          /* forms: 0 p, 1 p-q, 2 u, 3 q, 4 Q */
        static const int s8[5] = {1, 1, 2, 2, 2}, e8[5] = {4, 4, 4, 2, 1};
        nf = 5; for (int m = 0; m < 5; m++) { side[m] = s8[m]; ex[m] = e8[m]; }
    } else if (fam == FAM_Z12) {  /* forms: 0 q, 1 u, 2 A2, 3 A6, 4 p, 5 p-q, 6 A3 */
        static const int s12[7] = {2, 2, 2, 2, 1, 1, 1}, e12[7] = {2, 6, 3, 1, 6, 6, 2};
        nf = 7; for (int m = 0; m < 7; m++) { side[m] = s12[m]; ex[m] = e12[m]; }
    } else return -21;
    u64 B_[MAXB]; int cl[MAXB], vv[MAXB], n = 2;
    B_[0] = (u64)-1; B_[1] = 2; cl[0] = cl[1] = 0;
    for (int m = 0; m < nf; m++) for (int i = 0; i < ff[m].n; i++) {
        u64 l = ff[m].p[i]; if (l == 2) continue;
        int j = -1; for (int k = 2; k < n; k++) if (B_[k] == l) j = k;
        if (j < 0) { if (n >= MAXB) return -10; j = n++; B_[j] = l; cl[j] = 0; vv[j] = 0; }
        cl[j] |= side[m]; vv[j] += ex[m] * ff[m].e[i];
    }
    /* signs: B > 0 in both families; sign(B') = sign(Q) (z8) or sign(A6) (z12); sign(A) exactly in 128 bits */
    __int128 P = p, Q = q, A;
    if (fam == FAM_Z8) A = 8 * P * P * P * P - 16 * Q * P * P * P + 16 * Q * Q * P * P - 8 * Q * Q * Q * P + Q * Q * Q * Q;
    else { __int128 P2 = P * P, P4 = P2 * P2, Q2 = Q * Q, Q4 = Q2 * Q2;
        A = 24 * P4 * P4 - 96 * Q * P4 * P2 * P + 216 * Q2 * P4 * P2 - 312 * Q2 * Q * P4 * P + 288 * Q4 * P4 - 168 * Q4 * Q * P2 * P + 60 * Q4 * Q2 * P2 - 12 * Q4 * Q2 * Q * P + Q4 * Q4; }
    int sA = A > 0 ? 1 : (A < 0 ? -1 : 0);
    if (sA == 0) return -13;
    i64 Fs = (fam == FAM_Z8) ? (8 * p * p - 8 * p * q + q * q) : (6 * p * p - 6 * p * q + q * q);
    int sB = 1, sB2 = Fs > 0 ? 1 : -1;
    int neg1 = (sB < 0) || (sB2 > 0 && sA > 0);
    int neg2 = (sB2 < 0) || (sB > 0 && sA < 0);
    u64 R1[3 * MAXB + 8], R2[3 * MAXB + 8]; int n1 = 0, n2 = 0;
    u64 m1 = 1, m5 = 0;
    for (int k = 2; k < n; k++) { int r8 = (int)(B_[k] & 7); if ((r8 & 3) == 3) m1 |= 1ULL << k; if (r8 == 3 || r8 == 5) m5 |= 1ULL << k; }
    const u64 vb2 = 2;                                    /* the v-bit at 2 is carried only by the basis element 2 */
    if (!neg1) R1[n1++] = 1;
    if (!neg2) R2[n2++] = 1;
    /* l = 2 */
    {   int iq = (fam == FAM_Z8) ? 3 : 0, iu = (fam == FAM_Z8) ? 2 : 1;
        int v2q = 0, v2u = 0;
        for (int i = 0; i < ff[iq].n; i++) if (ff[iq].p[i] == 2) v2q = ff[iq].e[i];
        for (int i = 0; i < ff[iu].n; i++) if (ff[iu].p[i] == 2) v2u = ff[iu].e[i];
        if (v2q == 0) { R2[n2++] = vb2; R2[n2++] = m1; R2[n2++] = m5; }          /* W(E) all, W(E') 0 */
        else if (fam == FAM_Z12) {
            if (v2q == 1) { R1[n1++] = vb2; R1[n1++] = m1; R1[n1++] = m5; }       /* 0, all */
            else { R1[n1++] = vb2; R1[n1++] = m1; R2[n2++] = vb2; }               /* <5>, units */
        } else {
            if (v2q == 1 && v2u == 2) { R1[n1++] = vb2; R1[n1++] = m1; R2[n2++] = vb2; }   /* <5>, units */
            else if (v2q == 1) { R1[n1++] = vb2; R1[n1++] = m1; R1[n1++] = m5; }            /* 0, all */
            else if (v2q == 2) { R1[n1++] = vb2; R1[n1++] = m1 ^ m5; R2[n2++] = vb2 ^ m1; } /* <3>, <5,-2> */
            else { R1[n1++] = vb2; R1[n1++] = m5; R2[n2++] = m1; }                          /* <-1>, <5,2> */
        }
    }
    for (int j = 2; j < n; j++) {
        u64 l = B_[j];
        u64 lm = 0;
        if ((l & 3) == 3) lm |= 1;
        { int r8 = (int)(l & 7); if (r8 == 3 || r8 == 5) lm |= 2; }
        for (int k = 2; k < n; k++) if (k != j && nonres(B_[k], l, T)) lm |= 1ULL << k;
        u64 vrow = 1ULL << j;
        if (fam == FAM_Z12 && l == 3) {
            i64 p3 = ((p % 3) + 3) % 3, q3 = ((q % 3) + 3) % 3;
            if (q3 == 0) { R1[n1++] = vrow ^ lm; R2[n2++] = lm; }                 /* <3*2>, <3> */
            else if (p3 == 0 || p3 == q3) { R2[n2++] = vrow; R2[n2++] = lm; }   /* 3 | p(p-q): all, 0 */
            else { R1[n1++] = vrow; R1[n1++] = lm; }                           /* 3 | u: 0, all */
            continue;
        }
        if (cl[j] == 3) return -14;
        /* A mod l */
        u64 pm = (u64)((p % (i64)l + (i64)l) % (i64)l), qm = (u64)((q % (i64)l + (i64)l) % (i64)l), Am;
        #define MM(a, b) mulmod_any((a), (b), l)
        if (fam == FAM_Z8) { u64 p2 = MM(pm, pm), p3 = MM(p2, pm), p4 = MM(p3, pm), q2 = MM(qm, qm), q3 = MM(q2, qm), q4 = MM(q3, qm);
            Am = (MM(8, p4) + l - MM(16, MM(qm, p3)) + MM(16, MM(q2, p2)) + l - MM(8, MM(q3, pm)) + q4) % l; }
        else { u64 pw[9], qw[9]; pw[0] = qw[0] = 1; for (int k = 1; k <= 8; k++) { pw[k] = MM(pw[k - 1], pm); qw[k] = MM(qw[k - 1], qm); }
            static const i64 c12[9] = {24, -96, 216, -312, 288, -168, 60, -12, 1};    /* sum c_k q^k p^(8-k) */
            u64 acc = 0; for (int k = 0; k <= 8; k++) { u64 t = MM((u64)(c12[k] < 0 ? -c12[k] : c12[k]) % l, MM(qw[k], pw[8 - k])); acc = c12[k] < 0 ? (acc + l - t) % l : (acc + t) % l; }
            Am = acc; }
        #undef MM
        if (cl[j] == 1) {
            int split = legendre(Am, l) == 1;
            if (split || (vv[j] & 1)) { R2[n2++] = vrow; R2[n2++] = lm; }
            else { R1[n1++] = vrow; R2[n2++] = vrow; }
        } else {
            int split = legendre(mulmod_any(l - Am, 2, l), l) == 1;
            if (split || (vv[j] & 1)) { R1[n1++] = vrow; R1[n1++] = lm; }
            else { R1[n1++] = vrow; R2[n2++] = vrow; }
        }
    }
    int rk[2] = {0, 0};
    for (int g = 0; g < 2; g++) {
        u64 *R = g ? R2 : R1; int nr = g ? n2 : n1; u64 piv[64]; u64 pmk = 0;
        for (int i = 0; i < nr; i++) { u64 y = R[i];
            while (y) { int c = ctz64(y); if (!((pmk >> c) & 1)) { piv[c] = y; pmk |= 1ULL << c; rk[g]++; break; } y ^= piv[c]; } }
    }
    int s1 = n - rk[0], s2 = n - rk[1];
    if (sd) { sd[0] = s1; sd[1] = s2; }
    if (Bout) { for (int k = 0; k < n; k++) Bout[k] = B_[k]; *nout = n;
                for (int i = 0; i < n1; i++) R1out[i] = R1[i]; *n1out = n1; for (int i = 0; i < n2; i++) R2out[i] = R2[i]; *n2out = n2; }
    return s1 + s2 - 2;
}

/* bound for one pair: min over the family's curves (z8: E1 only unless all) of dim Sel2 - 2.
 * returns 0..250 bound, 255 singular fibre, 254 factorisation failure, 253 other error. */
/* target > 0: early exit allowed; the result is exact when >= target, and some value < target otherwise. */
/* all: bit 0 = z8 also E2;  bit 1 = add the 3-isogeny bound B3 (z2z6, z12): result min(B2, B3), dims[1] = B3 (or -1 if skipped);
   bit 2 = a parity is required (bit 3 = that parity): with bit 1 and target > 0, pairs of the wrong parity exit after factoring */
HDN int pair_bound(int fam, i64 p, i64 q, int all, const u32 *spf, u64 lim, const seltab_t *T, int *dims, int target) {
    i64 F[MAXF];
    {   /* domain: quadratic forms of z8, z12, z2z8 are at most 17 h^2 < 2^63 for h <= 7e8; linear z2z6: 10 h */
        i64 ap = p < 0 ? -p : p, aq = q < 0 ? -q : q, h = ap > aq ? ap : aq;
        if (h > (fam == FAM_Z2Z6 ? (i64)900000000000000000LL : fam == FAM_Z9 ? (i64)900000LL : (i64)700000000LL)) return SEL_ERROR;
    }
    int nf = family_forms(fam, p, q, F);
    if (nf <= 0) return SEL_SINGULAR;
    fact_t ff[MAXF];
    int b3first = (all & 2) && (fam == FAM_Z2Z6 || fam == FAM_Z12) && target > 0;
    int b3v = 1000;
    if (b3first) {                                /* all forms are used by these families: factor, then B3 before any curve setup */
        for (int m = 0; m < nf; m++) { u64 a = F[m] < 0 ? (u64)(-F[m]) : (u64)F[m]; if (factor64(a, &ff[m], spf, lim)) return SEL_FACTFAIL; }
        b3v = b3_fast(fam, p, q, ff, target, (all & 4) ? ((all >> 3) & 1) : -1);
        if (b3v < -1) return SEL_ERROR;
        if (b3v < target) return b3v;
    }
    if (fam == FAM_Z9) {
        for (int m = 0; m < nf; m++) { u64 a = F[m] < 0 ? (u64)(-F[m]) : (u64)F[m]; if (factor64(a, &ff[m], spf, lim)) return SEL_FACTFAIL; }
        int b0 = b3_z9(p, q, ff, 0, target); if (b0 < -1) return SEL_ERROR;
        if (dims) dims[0] = b0;
        if (target > 0 && b0 < target) return b0;
        int b1 = b3_z9(p, q, ff, 1, target); if (b1 < -1) return SEL_ERROR;
        if (dims) dims[1] = b1;
        return b0 < b1 ? b0 : b1;
    }
    if (fam == FAM_Z10) {
        for (int m = 0; m < nf; m++) { u64 a = F[m] < 0 ? (u64)(-F[m]) : (u64)F[m]; if (factor64(a, &ff[m], spf, lim)) return SEL_FACTFAIL; }
        int sd[2]; int b = bphi_z10(p, q, ff, T, sd);
        if (dims) { dims[0] = sd[0]; dims[1] = sd[1]; }
        return b < 0 ? SEL_ERROR : b;
    }
    curve_t cv[3];
    int nc = family_curves(fam, F, nf, cv);
    if (fam == FAM_Z8 && !(all & 1)) nc = 1;
    if (!b3first) for (int m = 0; m < nf; m++) {
        int used = 0;
        for (int c = 0; c < nc; c++) for (int t = 0; t < 3; t++) if (cv[c].e[t][m]) used = 1;
        if (!used) { ff[m].n = 0; continue; }
        u64 a = F[m] < 0 ? (u64)(-F[m]) : (u64)F[m];
        if (factor64(a, &ff[m], spf, lim)) return SEL_FACTFAIL;
    }
#ifdef ONLYFACT
    return ff[0].n;
#endif
    if (fam == FAM_Z8 && target > 0 && (all & 4) && T->tabw2) {     /* parity cut from the factorisation */
        int pr = z8_parity(&cv[0], ff, T->tabw2);
        if (pr >= 0 && pr != ((all >> 3) & 1)) return target - 1;
    }
    int best = b3first ? b3v : 1000;
    for (int c = 0; c < nc; c++) {
        int d = sel2_dim(&cv[c], ff, T, target > 0 ? target + 2 : 0);
        if (dims) dims[c] = d;
        if (d < 0) return SEL_ERROR;
        if (d - 2 < best) best = d - 2;
        if (target > 0 && best < target) return best;
    }
    if ((all & 2) && (fam == FAM_Z2Z6 || fam == FAM_Z12) && !b3first) {
        if (dims) dims[1] = -1;
        if (target <= 0 || best >= target) {
            int b3 = b3_fam(fam, p, q, ff);
            if (b3 < -1) return SEL_ERROR;
            if (dims) dims[1] = b3;
            if (b3 < best) best = b3;
        }
    }
    return best;
}

/* host-side construction of the QR bitmaps (call once; arrays sized by qr_sizes) */
static inline void qr_build(u32 *qroff, u64 *qrbits, u64 nwords) {
    u64 off = 0;
    for (u32 b = 0; b < QRMAX; b++) qroff[b] = 0;
    for (u64 i = 0; i < nwords; i++) qrbits[i] = 0;
    for (u32 b = 3; b < QRMAX; b += 2) {
        int pr = 1; for (u32 d = 3; d * d <= b; d += 2) if (b % d == 0) { pr = 0; break; }
        if (!pr) continue;
        qroff[b] = (u32)off;
        for (u64 x = 1; x < b; x++) { u64 r = x * x % b, o = off + r; qrbits[o >> 6] |= 1ULL << (o & 63); }
        off += b;
    }
}
#define QR_WORDS ((u64)(QRMAX * QRMAX / 2 / 64 + 64))
