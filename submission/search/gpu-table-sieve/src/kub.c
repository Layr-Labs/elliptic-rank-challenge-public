// kub.c - Mestre-Nagao table sieve for the Tate-normal-form torsion families Z/5, Z/6, Z/2xZ/4
// (own code; method: Mestre 1982 / Nagao 1992 sums S(X) = sum_{l<=X} -a_l log(l)/l, per-prime tables indexed by t mod l)
//
// build:  cc -O3 -march=native -o kub kub.c -lm -lpthread
// usage:  kub sieve FAM Hmax qlo qhi thr1 thr2 thr3 [threads]     (FAM = 5 | 6 | 42)
//         kub check FAM p q                                      (print S300/S2000/S10000 and H for one parameter t=p/q)
// FAM 5:  y^2+(1-t)xy-t y = x^3-t x^2                      (Z/5)
// FAM 6:  y^2+(1-t)xy-(t+t^2)y = x^3-(t+t^2)x^2             (Z/6)
// FAM 42: y^2+xy-b y = x^3-b x^2, b = v^2-1/16, v=t         (Z/2xZ/4)
// FAM 44: y^2 = x (x+p^2) (x+q^2), (p,q) = (b,a) with a > b >= 1 coprime  (Z/2xZ/4 in Legendre form: no duplicates, no 2-power scaling overhead; a4 = p^2 q^2 != 0)
// sieve output (stdout): p q S300 S2000 S10000 H   for every t=p/q (gcd=1, q in [qlo,qhi]) with H<Hmax and S300>=thr1,
//                        S2000>=thr2, S10000>=thr3.  progress / counters go to stderr.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <time.h>
static double now(void){struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+1e-9*t.tv_nsec;}
#if defined(__ARM_NEON) && !defined(__CUDACC__)
#include <arm_neon.h>
#endif

#define PMAX 10000
#define S1PMAX 300
#define S2PMAX 2000
#define SCALE 8.0
#define BLK 512

static int FAM;
static int SYM = 1;   // restrict to the fundamental domain of the t -> -1/t (Z/5) and v -> -v (Z/2xZ/4) symmetries
static int nprimes, n1, n2;              // odd primes < PMAX, < S1PMAX, < S2PMAX
static int primes[1300];
static int8_t *tab[1300];                // tab[i][tau], tau in [0,l)
static uint16_t *inv[1300];              // inv[i][r] = r^-1 mod l (0 if r==0)
static int8_t *pat[1300];                // pat[i][r*l + j] = tab[i][(j*inv[r]) % l], r in [1,l); (stage-1 primes only)

// ---------------------------------------------------------------- family models (table-driven, from fams.txt written by famgen.gp)
// a_i'(p,q) = q^(w_i*FD) * A_i(p/q),  [a1,a2,a3,a4,a6] with weights w = 1,2,3,4,6;  A_i integer polynomials in tau = p/q.
// fams.txt blocks "A1 .. A2 .. A3 .." are mandatory (a4 = a6 = 0 unless given); optional trailing blocks "A4 n c0 .." and "A6 n c0 .." (the old 3-block parsers skip them).
#define MAXC 16
#define NCO 5
static const int AWT[NCO] = {1, 2, 3, 4, 6};
static int FD;                       // scale degree d
static int NAC[NCO];
static long long AC[NCO][MAXC];      // ascending coefficients
static void load_family(int fam) {
    const char *fn = getenv("KUB_FAMS"); if (!fn) fn = "fams.txt";
    FILE *f = fopen(fn, "r");
    if (!f) { fprintf(stderr, "cannot open family table %s\n", fn); exit(3); }
    char tok[64]; int found = 0;
    for (int i = 0; i < NCO; i++) { NAC[i] = 1; AC[i][0] = 0; }
    while (fscanf(f, "%63s", tok) == 1) {
        if (strcmp(tok, "FAM")) continue;
        char name[64]; int d;
        if (fscanf(f, "%63s %d", name, &d) != 2) break;
        if (atoi(name) != fam) continue;
        FD = d;
        for (int i = 0; i < 3; i++) {
            char an[16]; int n;
            if (fscanf(f, "%15s %d", an, &n) != 2 || n > MAXC) { fprintf(stderr, "bad family table\n"); exit(3); }
            NAC[i] = n;
            for (int j = 0; j < n; j++) if (fscanf(f, "%lld", &AC[i][j]) != 1) { fprintf(stderr, "bad family table\n"); exit(3); }
        }
        for (int extra = 0; extra < 2; extra++) {            // optional A4 / A6 blocks
            long pos = ftell(f); char an[16]; int n;
            if (fscanf(f, "%15s", an) != 1) break;
            int slot = !strcmp(an, "A4") ? 3 : (!strcmp(an, "A6") ? 4 : -1);
            if (slot < 0) { fseek(f, pos, SEEK_SET); break; }
            if (fscanf(f, "%d", &n) != 1 || n > MAXC) { fprintf(stderr, "bad family table\n"); exit(3); }
            NAC[slot] = n;
            for (int j = 0; j < n; j++) if (fscanf(f, "%lld", &AC[slot][j]) != 1) { fprintf(stderr, "bad family table\n"); exit(3); }
        }
        found = 1; break;
    }
    fclose(f);
    if (!found) { fprintf(stderr, "family %d not in %s\n", fam, fn); exit(3); }
}
// [a1,a2,a3,a4,a6](tau) mod l for tau in [0,l)
static void coeffs_mod(int fam, int64_t tau, int64_t l, int64_t *a) {
    (void)fam;
    for (int i = 0; i < NCO; i++) {
        int64_t acc = 0;
        for (int j = NAC[i] - 1; j >= 0; j--) {
            int64_t c = AC[i][j] % l; if (c < 0) c += l;
            acc = (acc * tau + c) % l;
        }
        a[i] = acc;
    }
}

// character sum  sum_x chi(4x^3 + b2 x^2 + 2 b4 x + b6) over F_l, chi = Legendre symbol (table chi[0..l));  b2 = a1^2+4a2, b4 = a1 a3 + 2 a4, b6 = a3^2 + 4 a6
static int charsum(const int8_t *chi, int64_t l, const int64_t *a) {
    const int64_t a1 = a[0], a2 = a[1], a3 = a[2], a4 = a[3], a6 = a[4];
    int64_t b2 = (a1 * a1 + 4 * a2) % l, b4 = (a1 * a3 + 2 * a4) % l, b6 = (a3 * a3 + 4 * a6) % l;
    int64_t c3 = 4 % l, c2 = b2, c1 = (2 * b4) % l, c0 = b6;
    int s = 0;
    for (int64_t x = 0; x < l; x++) {
        int64_t f = ((c3 * x + c2) % l * x + c1) % l * x + c0;
        s += chi[f % l];
    }
    return s;
}

typedef struct { int lo, hi; } rng_t;
static void *build_tables_thr(void *arg) {
    rng_t *r = (rng_t *)arg;
    for (int i = r->lo; i < r->hi; i++) {
        int64_t l = primes[i];
        int8_t *chi = (int8_t *)calloc(l, 1);
        for (int64_t x = 1; x < l; x++) chi[x * x % l] = 1;                // mark squares
        for (int64_t x = 1; x < l; x++) chi[x] = chi[x] ? 1 : -1;          // chi(x)
        chi[0] = 0;
        tab[i] = (int8_t *)malloc(l);
        double lw = log((double)l) / (double)l;
        for (int64_t tau = 0; tau < l; tau++) {
            int64_t ac[NCO];
            coeffs_mod(FAM, tau, l, ac);
            int cs = charsum(chi, l, ac);                                 // = -a_l
            double w = cs * lw * SCALE;
            int v = (int)floor(w + 0.5);
            if (v > 127) v = 127;
            if (v < -127) v = -127;
            tab[i][tau] = (int8_t)v;
        }
        free(chi);
        inv[i] = (uint16_t *)calloc(l, sizeof(uint16_t));
        for (int64_t a = 1; a < l; a++) { // inverse by Fermat
            int64_t b = 1, base = a, e = l - 2;
            while (e) { if (e & 1) b = b * base % l; base = base * base % l; e >>= 1; }
            inv[i][a] = (uint16_t)b;
        }
    }
    return NULL;
}

static void init_primes(void) {
    static char comp[PMAX + 1];
    nprimes = 0;
    for (int i = 3; i < PMAX; i += 2) {
        if (comp[i]) continue;
        primes[nprimes++] = i;
        for (int j = i * i; j < PMAX; j += 2 * i) comp[j] = 1;
    }
    n1 = n2 = 0;
    for (int i = 0; i < nprimes; i++) { if (primes[i] < S1PMAX) n1 = i + 1; if (primes[i] < S2PMAX) n2 = i + 1; }
}

static int load_tab_cache(void) {
    char fn[512];
    const char *dir = getenv("KUB_TABDIR");
    if (!dir) return 0;
    snprintf(fn, sizeof fn, "%s/tab_%d.bin", dir, FAM);
    FILE *f = fopen(fn, "rb");
    if (!f) return 0;
    for (int i = 0; i < nprimes; i++) {
        tab[i] = (int8_t *)malloc(primes[i]);
        if (fread(tab[i], 1, primes[i], f) != (size_t)primes[i]) { fclose(f); return 0; }
    }
    fclose(f);
    return 1;
}
static void save_tab_cache(void) {
    char fn[512];
    const char *dir = getenv("KUB_TABDIR");
    if (!dir) return;
    snprintf(fn, sizeof fn, "%s/tab_%d.bin", dir, FAM);
    FILE *f = fopen(fn, "wb");
    if (!f) return;
    for (int i = 0; i < nprimes; i++) fwrite(tab[i], 1, primes[i], f);
    fclose(f);
}

static void build_tables(int nthr) {
    init_primes();
    if (load_tab_cache()) {
        for (int i = 0; i < nprimes; i++) {
            int64_t l = primes[i];
            inv[i] = (uint16_t *)calloc(l, sizeof(uint16_t));
            for (int64_t a = 1; a < l; a++) {
                int64_t b = 1, base = a, e = l - 2;
                while (e) { if (e & 1) b = b * base % l; base = base * base % l; e >>= 1; }
                inv[i][a] = (uint16_t)b;
            }
        }
    } else {
        pthread_t th[64]; rng_t rg[64];
        if (nthr > 64) nthr = 64;
        double tot = 0; for (int i = 0; i < nprimes; i++) tot += (double)primes[i] * primes[i];
        int idx = 0;
        for (int t = 0; t < nthr; t++) {
            rg[t].lo = idx; double acc = 0;
            while (idx < nprimes && (acc < tot / nthr || t == nthr - 1)) { acc += (double)primes[idx] * primes[idx]; idx++; }
            rg[t].hi = idx;
            pthread_create(&th[t], NULL, build_tables_thr, &rg[t]);
        }
        for (int t = 0; t < nthr; t++) pthread_join(th[t], NULL);
        save_tab_cache();
    }
    for (int i = 0; i < n1; i++) {
        int64_t l = primes[i];
        pat[i] = (int8_t *)malloc(l * l);
        for (int64_t r = 1; r < l; r++)
            for (int64_t j = 0; j < l; j++) pat[i][r * l + j] = tab[i][(j * inv[i][r]) % l];
        for (int64_t j = 0; j < l; j++) pat[i][j] = 0;                    // r == 0 (l | q): no contribution
    }
}

// ---------------------------------------------------------------- height (approximate, from the integral model)
static void model(int fam, double p, double q, double *a) {      // a[0..4] = a1,a2,a3,a4,a6 of the integral model at (p,q)
    (void)fam;
    for (int i = 0; i < NCO; i++) {
        double acc = 0, pw = 1;                 // sum_j c_j p^j q^(w-j), w = AWT[i]*FD
        int w = AWT[i] * FD;
        double qp[64]; qp[0] = 1; for (int k = 1; k <= w; k++) qp[k] = qp[k-1] * q;
        for (int j = 0; j < NAC[i]; j++) { if (w - j >= 0) acc += (double)AC[i][j] * pw * qp[w - j]; pw *= p; }
        a[i] = acc;
    }
}
static void bvals(const double *a, double *b2, double *b4, double *b6) { *b2 = a[0] * a[0] + 4 * a[1]; *b4 = a[0] * a[2] + 2 * a[3]; *b6 = a[2] * a[2] + 4 * a[4]; }
static double hnaive(int fam, double p, double q) {
    double a[NCO]; model(fam, p, q, a);
    double b2, b4, b6; bvals(a, &b2, &b4, &b6);
    double c4 = b2 * b2 - 24 * b4, c6 = -b2 * b2 * b2 + 36 * b2 * b4 - 216 * b6;
    double x = fabs(c4), y = fabs(c6);
    double l4 = x > 0 ? 3 * log(x) : -1e300, l6 = y > 0 ? 2 * log(y) : -1e300;
    return (l4 > l6 ? l4 : l6) / 12.0;
}

// S sums for one (p,q) by scalar code; staged with early exit (thresholds t2,t3; pass t2=t3=-1e9 for no exit)
static int EARLY = 1;       // 1: early exits at S(600), S(1000), S(4000) (CPU sieve only, wide margins); the GPU driver sets 0 = exact thresholds only
static int sums_staged(int64_t p, int64_t q, double t2, double t3, double *s300, double *s2000, double *s10000) {
    int a = 0, b = 0, c = 0;
    int i1000 = 0, i600 = 0, i4000 = 0;
    for (int i = 0; i < nprimes; i++) { if (primes[i] < 600) i600 = i + 1; if (primes[i] < 1000) i1000 = i + 1; if (primes[i] < 4000) i4000 = i + 1; }
    for (int i = 0; i < nprimes; i++) {
        int64_t l = primes[i];
        int64_t qm = q % l;
        if (qm != 0) {
            int64_t pm = ((p % l) + l) % l;
            int v = tab[i][(pm * inv[i][qm]) % l];
            c += v; if (i < n2) b += v; if (i < n1) a += v;
        }
        // The old margins (t2-6 at 600, t2-4 at 1000, t3-8 at 4000) dropped about a third of the frozen rank-8/9 records (S600 = 19.75, 17.0, 19.75 < 20).
        // Now: wide margins (the profile of the rank-7..9 curves has S2000 - S600 >= 4.4 and S10000 - S4000 >= 6) and only when EARLY is set.
        if (EARLY && i + 1 == i600 && b / SCALE < t2 - 14.0) return 0;
        if (EARLY && i + 1 == i1000 && b / SCALE < t2 - 12.0) return 0;
        if (i + 1 == n2) { if (b / SCALE < t2) return 0; }
        if (EARLY && i + 1 == i4000 && c / SCALE < t3 - 14.0) return 0;
    }
    *s300 = a / SCALE; *s2000 = b / SCALE; *s10000 = c / SCALE;
    return 1;
}
static void sums_scalar(int64_t p, int64_t q, double *s300, double *s2000, double *s10000) {
    sums_staged(p, q, -1e9, -1e9, s300, s2000, s10000);
}

// ---------------------------------------------------------------- sieve
static double Hmax, Hmin = -1e9, thr1, thr2, thr3;
static int64_t qlo, qhi, qnext;
static pthread_mutex_t qmu = PTHREAD_MUTEX_INITIALIZER, omu = PTHREAD_MUTEX_INITIALIZER;
static volatile uint64_t cnt_cand = 0, cnt_s1 = 0, cnt_s2 = 0, cnt_s3 = 0;

static int64_t gcd64(int64_t a, int64_t b) { if (a < 0) a = -a; while (b) { int64_t t = a % b; a = b; b = t; } return a; }

// Phi(tau) = max(|c4|^3, c6^2) of the homogeneous model at (p,q)=(tau,1); H(p,q) = ln q + ln Phi(p/q) / 12
static double phi_of(double tau) {
    double a[NCO]; model(FAM, tau, 1.0, a);
    double b2, b4, b6; bvals(a, &b2, &b4, &b6);
    double c4 = b2 * b2 - 24 * b4, c6 = -b2 * b2 * b2 + 36 * b2 * b4 - 216 * b6;
    double x = fabs(c4), y = fabs(c6);
    double u = x * x * x, v = y * y;
    return u > v ? u : v;
}
#define NGRID 6200
static double gx[NGRID], gphi[NGRID];
static int ngrid;
static void init_grid(void) {
    // ascending: -1e10 ... -1e-3, 0, 1e-3 ... 1e10 geometric ratio 1.01
    int n = 0;
    double r = 1.01;
    int m = 0; for (double x = 1e-3; x < 1e10; x *= r) m++;
    double *pos = (double *)malloc(m * sizeof(double)); m = 0;
    for (double x = 1e-3; x < 1e10; x *= r) pos[m++] = x;
    for (int i = m - 1; i >= 0; i--) gx[n++] = -pos[i];
    gx[n++] = 0.0;
    for (int i = 0; i < m; i++) gx[n++] = pos[i];
    ngrid = n; free(pos);
    for (int i = 0; i < n; i++) gphi[i] = phi_of(gx[i]);
}
// largest q that can still have candidates with H < Hlim:  q <= exp(Hlim) / Phi_min^(1/12)  (+5 % margin)
static int64_t q_auto_cap(double Hlim) {
    double pm = 1e300;
    for (int i = 0; i < ngrid; i++) if (gphi[i] < pm) pm = gphi[i];
    return (int64_t)(pow(exp(Hlim) * pow(pm, -1.0 / 12.0), 1.0 / FD) * 1.05) + 2;
}

// intervals of p (inclusive) with H(p,q) < Hlim (plus a small slack); returns count (<= 8)
static int p_intervals(int64_t q, double Hlim, int64_t *lo, int64_t *hi) {
    double L = exp(12.0 * (Hlim - FD * log((double)q)));
    int cnt = 0;
    int i = 0;
    while (i < ngrid && cnt < 8) {
        while (i < ngrid && gphi[i] >= L) i++;
        if (i >= ngrid) break;
        int j = i;
        while (j + 1 < ngrid && gphi[j + 1] < L) j++;
        double t1 = gx[i], t2 = gx[j];
        if (i > 0) { // refine left edge between gx[i-1] (out) and gx[i] (in)
            double a = gx[i - 1], b = gx[i];
            for (int it = 0; it < 50; it++) { double m = 0.5 * (a + b); if (phi_of(m) < L) b = m; else a = m; }
            t1 = a;
        }
        if (j + 1 < ngrid) {
            double a = gx[j], b = gx[j + 1];
            for (int it = 0; it < 50; it++) { double m = 0.5 * (a + b); if (phi_of(m) < L) a = m; else b = m; }
            t2 = b;
        }
        lo[cnt] = (int64_t)floor(t1 * (double)q) - 1;
        hi[cnt] = (int64_t)ceil(t2 * (double)q) + 1;
        cnt++;
        i = j + 1;
    }
    return cnt;
}

// intervals of p with Hlo <= H(p,q) < Hhi  (= intervals(Hhi) minus intervals(Hlo)); returns count (<= 16)
static int p_shell(int64_t q, double Hlo, double Hhi, int64_t *lo, int64_t *hi) {
    int64_t alo[8], ahi[8], blo[8], bhi[8];
    int na = p_intervals(q, Hhi, alo, ahi);
    int nb = Hlo > -1e8 ? p_intervals(q, Hlo, blo, bhi) : 0;
    for (int j = 0; j < nb; j++) { blo[j] += 3; bhi[j] -= 3; }   // shrink: the interval slack must stay in the shell (host re-checks H >= Hmin)
    int n = 0;
    for (int i = 0; i < na; i++) {
        int64_t cur = alo[i];
        for (int j = 0; j < nb; j++) {
            if (blo[j] > bhi[j] || bhi[j] < cur || blo[j] > ahi[i]) continue;
            if (blo[j] > cur && n < 16) { lo[n] = cur; hi[n] = blo[j] - 1; n++; }
            if (bhi[j] + 1 > cur) cur = bhi[j] + 1;
        }
        if (cur <= ahi[i] && n < 16) { lo[n] = cur; hi[n] = ahi[i]; n++; }
    }
    return n;
}

typedef struct { int tid; } targ_t;

static void *sieve_thr(void *arg) {
    (void)arg;
    int8_t *ext[1300];
    for (int i = 0; i < n1; i++) ext[i] = (int8_t *)aligned_alloc(64, ((BLK + primes[i] + 64) + 63) / 64 * 64);
    int8_t *S = (int8_t *)aligned_alloc(64, BLK + 64);
    int8_t thr1i = (int8_t)floor(thr1 * SCALE + 0.5);
    for (;;) {
        pthread_mutex_lock(&qmu);
        int64_t q = qnext++;
        pthread_mutex_unlock(&qmu);
        if (q > qhi) break;
        int64_t ilo[16], ihi[16];
        int nint = p_shell(q, Hmin, Hmax, ilo, ihi);
        if (nint == 0) continue;
        if (SYM) { // clip to the fundamental domain
            int m = 0;
            for (int iv = 0; iv < nint; iv++) {
                int64_t lo = ilo[iv], hi = ihi[iv];
                if (FAM == 5) { if (lo < -q) lo = -q; if (hi > q) hi = q; }
                else if (FAM == 42) { if (lo < 1) lo = 1; }
                else if (FAM == 8) { if (hi > (q >> 1)) hi = q >> 1; }          // Z/8: t and 1-t give the same curve -> fundamental domain p <= q/2
                else if (FAM == 44) { if (lo < 1) lo = 1; if (hi > q - 1) hi = q - 1; }   // Legendre Z/2xZ/4: (a,b) ~ (b,a) ~ (-a,b): fundamental domain 1 <= p=b < q=a
                if (lo <= hi) { ilo[m] = lo; ihi[m] = hi; m++; }
            }
            nint = m;
            if (nint == 0) continue;
        }
        // build extended pattern rows for this q (rows replicated to length BLK + l)
        for (int i = 0; i < n1; i++) {
            int64_t l = primes[i], r = q % l;
            const int8_t *row = pat[i] + r * l;
            int64_t len = BLK + l, k = 0;
            while (k < len) { int64_t c = len - k < l ? len - k : l; memcpy(ext[i] + k, row, c); k += c; }
        }
        uint64_t local_c = 0, local_s1 = 0, local_s2 = 0, local_s3 = 0;
        char *obuf = NULL; size_t olen = 0, ocap = 0;
        for (int iv = 0; iv < nint; iv++)
        for (int64_t p0 = ilo[iv]; p0 <= ihi[iv]; p0 += BLK) {
            int64_t len = ihi[iv] - p0 + 1; if (len > BLK) len = BLK;
            {
                int offs[1300];
                for (int i = 0; i < n1; i++) { int64_t l = primes[i]; offs[i] = (int)(((p0 % l) + l) % l); }
#if defined(__ARM_NEON) && !defined(__CUDACC__)
                for (int tile = 0; tile < BLK; tile += 128) {
                    int8x16_t s0 = vdupq_n_s8(0), s1 = s0, s2 = s0, s3 = s0, s4 = s0, s5 = s0, s6 = s0, s7 = s0;
                    for (int i = 0; i < n1; i++) {
                        const int8_t *e = ext[i] + offs[i] + tile;
                        s0 = vqaddq_s8(s0, vld1q_s8(e));       s1 = vqaddq_s8(s1, vld1q_s8(e + 16));
                        s2 = vqaddq_s8(s2, vld1q_s8(e + 32));  s3 = vqaddq_s8(s3, vld1q_s8(e + 48));
                        s4 = vqaddq_s8(s4, vld1q_s8(e + 64));  s5 = vqaddq_s8(s5, vld1q_s8(e + 80));
                        s6 = vqaddq_s8(s6, vld1q_s8(e + 96));  s7 = vqaddq_s8(s7, vld1q_s8(e + 112));
                    }
                    vst1q_s8(S + tile, s0);       vst1q_s8(S + tile + 16, s1);  vst1q_s8(S + tile + 32, s2);  vst1q_s8(S + tile + 48, s3);
                    vst1q_s8(S + tile + 64, s4);  vst1q_s8(S + tile + 80, s5);  vst1q_s8(S + tile + 96, s6);  vst1q_s8(S + tile + 112, s7);
                }
#else
                memset(S, 0, BLK);
                for (int i = 0; i < n1; i++) {
                    const int8_t *e = ext[i] + offs[i];
                    for (int k = 0; k < BLK; k++) { int v = S[k] + e[k]; S[k] = v > 127 ? 127 : (v < -127 ? -127 : v); }
                }
#endif
            }
            local_c += (uint64_t)len;
            for (int k = 0; k < len; k++) {
                if (S[k] < thr1i) continue;
                int64_t p = p0 + k;
                if (p == 0 || gcd64(p, q) != 1) continue;
                double h = hnaive(FAM, (double)p, (double)q);
                if (h >= Hmax || h < Hmin) continue;
                local_s1++;
                double s3, s2, s1;
                if (!sums_staged(p, q, thr2, thr3, &s1, &s2, &s3)) continue;
                if (s1 < thr1 || s2 < thr2) continue;
                local_s2++;
                if (s3 < thr3) continue;
                local_s3++;
                char line[160];
                int n = snprintf(line, sizeof line, "%lld %lld %.2f %.2f %.2f %.4f\n", (long long)p, (long long)q, s1, s2, s3, h);
                if (olen + n + 1 > ocap) { ocap = ocap ? ocap * 2 : 4096; obuf = (char *)realloc(obuf, ocap); }
                memcpy(obuf + olen, line, n); olen += n;
            }
        }
        pthread_mutex_lock(&omu);
        if (olen) { fwrite(obuf, 1, olen, stdout); fflush(stdout); }
        cnt_cand += local_c; cnt_s1 += local_s1; cnt_s2 += local_s2; cnt_s3 += local_s3;
        pthread_mutex_unlock(&omu);
        free(obuf);
    }
    return NULL;
}

#ifndef KUB_NO_MAIN
int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: kub sieve FAM Hmax qlo qhi thr1 thr2 thr3 [threads] | kub check FAM p q\n"); return 1; }
    FAM = atoi(argv[2]);
    load_family(FAM);
    if (getenv("KUB_NOSYM")) SYM = 0;
    if (getenv("KUB_HMIN")) Hmin = atof(getenv("KUB_HMIN"));
    init_grid();
    if (!strcmp(argv[1], "check")) {
        build_tables(8);
        int64_t p = atoll(argv[3]), q = atoll(argv[4]);
        double a, b, c; sums_scalar(p, q, &a, &b, &c);
        printf("t=%lld/%lld  S300=%.2f S2000=%.2f S10000=%.2f  H=%.4f\n", (long long)p, (long long)q, a, b, c, hnaive(FAM, (double)p, (double)q));
        return 0;
    }
    if (!strcmp(argv[1], "profile")) {      // kub profile FAM FILE : for every line "p q ..." of FILE print the partial sums S(X) for X = 100 ... 10000
        build_tables(8);
        static const int cut[] = {100, 200, 300, 450, 600, 800, 1000, 1400, 2000, 4000, 10000};
        FILE *f = fopen(argv[3], "r"); if (!f) { perror("open"); return 1; }
        char line[512];
        while (fgets(line, sizeof line, f)) {
            long long p, q; if (sscanf(line, "%lld %lld", &p, &q) != 2) continue;
            int acc = 0, c = 0; printf("%lld %lld", p, q);
            for (int i = 0; i < nprimes && c < 11; i++) {
                int64_t l = primes[i]; int64_t qm = q % l;
                if (qm != 0) { int64_t pm = (((int64_t)p % l) + l) % l; acc += tab[i][(pm * inv[i][qm]) % l]; }
                while (c < 11 && (i + 1 == nprimes || primes[i + 1] > cut[c])) { printf(" %.2f", acc / SCALE); c++; }
            }
            printf("\n");
        }
        return 0;
    }
    Hmax = atof(argv[3]); qlo = atoll(argv[4]); qhi = atoll(argv[5]);
    thr1 = atof(argv[6]); thr2 = atof(argv[7]); thr3 = atof(argv[8]);
    { int64_t qa = q_auto_cap(Hmax); if (qhi > qa) qhi = qa; fprintf(stderr, "q range %lld..%lld (automatic cap %lld)\n", (long long)qlo, (long long)qhi, (long long)qa); }
    int nthr = argc > 9 ? atoi(argv[9]) : 8;
    double t0 = now();
    build_tables(nthr);
    double t1 = now();
    fprintf(stderr, "tables ready: %d odd primes <%d (stage1 %d, stage2 %d)\n", nprimes, PMAX, n1, n2);
    qnext = qlo;
    pthread_t th[64]; if (nthr > 64) nthr = 64;
    for (int t = 0; t < nthr; t++) pthread_create(&th[t], NULL, sieve_thr, NULL);
    for (int t = 0; t < nthr; t++) pthread_join(th[t], NULL);
    fprintf(stderr, "tables %.1fs  sieve %.1fs  rate %.3g candidates/s\n", t1 - t0, now() - t1, (double)cnt_cand / (now() - t1));
    fprintf(stderr, "done: candidates %llu  after S300 %llu  after S2000 %llu  output %llu\n",
            (unsigned long long)cnt_cand, (unsigned long long)cnt_s1, (unsigned long long)cnt_s2, (unsigned long long)cnt_s3);
    return 0;
}
#endif
