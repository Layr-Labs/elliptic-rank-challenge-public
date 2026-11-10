// tabx.cu - extended per-prime tables for the Mestre-Nagao sums beyond 10^4 (own code).
//
// For a Kubert family FAM (kub.c / fams.txt) and every prime l with PLO < l < PHI (default 10^4 .. 10^5) the table holds the integer a_l(E_tau) for all tau in F_l
// (tau = t mod l, E_tau = reduction of the integral model [a1,a2,a3,0,0] = [A1(tau), A2(tau), A3(tau), 0, 0]; a_l = -sum_x chi(4x^3 + b2 x^2 + 2 b4 x + b6)).
// With it the host can evaluate dS = S(10^5) - S(10^4) for a candidate in ~1 ms (one modular inverse and one lookup per prime) instead of 70 ms with PARI's ellap.
//
// build the tables (GPU, ~1-2 min):   nvcc -O3 -arch=sm_121 -DTABX_MAIN -o tabx_build tabx.cu       then   KUB_FAMS=fams.txt ./tabx_build build FAM tabdir
// check against a list "p q" :        ./tabx_build check FAM tabdir listfile        (prints the increment S(10^4 < l < 10^5) per line; compare with PARI)
// use in the sieve: kubg.cu includes this file with TABX_NO_MAIN and calls tabx_load(FAM) / tabx_sum(p, q).
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#ifndef TABX_NO_MAIN
#include <cuda_runtime.h>
#define KUB_NO_MAIN
#include "kub.c"
#endif

#define TABX_MAGIC 0x58424154u      // "TABX"
// file: uint32 magic, int32 fam, int32 plo, int32 phi, int32 nprimes, then per prime: uint32 l, int16 a[l]

// ---------------------------------------------------------------- host reader
static int tabx_n = 0;
static uint32_t *tabx_pr = NULL;
static const int16_t **tabx_tab = NULL;
static double *tabx_w = NULL;
static int tabx_ok = 0;

static int tabx_load(int fam) {
    const char *dir = getenv("KUB_TABDIR");
    if (!dir) return 0;
    char fn[512]; snprintf(fn, sizeof fn, "%s/tabx_%d.bin", dir, fam);
    int fd = open(fn, O_RDONLY);
    if (fd < 0) return 0;
    struct stat st; fstat(fd, &st);
    uint8_t *base = (uint8_t *)mmap(NULL, st.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
    close(fd);
    if (base == MAP_FAILED) return 0;
    uint32_t *h = (uint32_t *)base;
    if (h[0] != TABX_MAGIC || (int)h[1] != fam) { fprintf(stderr, "bad tabx file %s\n", fn); return 0; }
    tabx_n = (int)h[4];
    tabx_pr = (uint32_t *)malloc(sizeof(uint32_t) * tabx_n);
    tabx_tab = (const int16_t **)malloc(sizeof(int16_t *) * tabx_n);
    tabx_w = (double *)malloc(sizeof(double) * tabx_n);
    uint8_t *p = base + 5 * sizeof(uint32_t);
    for (int i = 0; i < tabx_n; i++) {
        uint32_t l = *(uint32_t *)p; p += sizeof(uint32_t);
        tabx_pr[i] = l; tabx_tab[i] = (const int16_t *)p; p += sizeof(int16_t) * (size_t)l;
        tabx_w[i] = log((double)l) / (double)l;
    }
    tabx_ok = 1;
    return 1;
}

static inline uint32_t tabx_modinv(uint32_t a, uint32_t m) {      // a in [1,m), m prime
    int64_t t0 = 0, t1 = 1; int64_t r0 = m, r1 = a;
    while (r1) { int64_t qq = r0 / r1; int64_t t2 = t0 - qq * t1; t0 = t1; t1 = t2; int64_t r2 = r0 - qq * r1; r0 = r1; r1 = r2; }
    return (uint32_t)(t0 < 0 ? t0 + m : t0);
}

// sum over the tabulated primes of -a_l log(l)/l for the parameter t = p/q  (primes dividing q contribute nothing, as in kub.c)
static double tabx_sum(int64_t p, int64_t q) {
    double s = 0;
    for (int i = 0; i < tabx_n; i++) {
        uint32_t l = tabx_pr[i];
        int64_t qm = q % (int64_t)l; if (qm < 0) qm += l;
        if (qm == 0) continue;
        int64_t pm = p % (int64_t)l; if (pm < 0) pm += l;
        uint32_t tau = (uint32_t)(((uint64_t)pm * tabx_modinv((uint32_t)qm, l)) % l);
        s -= (double)tabx_tab[i][tau] * tabx_w[i];
    }
    return s;
}

#ifdef TABX_MAIN
// ---------------------------------------------------------------- GPU builder
#define CK(x) do { cudaError_t e = (x); if (e != cudaSuccess) { fprintf(stderr, "CUDA error %s at %s:%d\n", cudaGetErrorString(e), __FILE__, __LINE__); exit(2); } } while (0)
__constant__ long long c_AC[NCO][MAXC];
__constant__ int c_NAC[NCO];

__global__ void __launch_bounds__(256) build_prime(int l, int16_t *out) {
    extern __shared__ uint32_t qr[];                         // bit table of the nonzero squares mod l
    const int nw = (l + 31) / 32;
    for (int i = threadIdx.x; i < nw; i += blockDim.x) qr[i] = 0;
    __syncthreads();
    for (int x = 1 + threadIdx.x; x <= (l - 1) / 2; x += blockDim.x) {
        unsigned long long sq = (unsigned long long)x * x % l;
        atomicOr(&qr[sq >> 5], 1u << (sq & 31));
    }
    __syncthreads();
    const int tau = blockIdx.x * blockDim.x + threadIdx.x;
    if (tau >= l) return;
    long long a[NCO];
    for (int i = 0; i < NCO; i++) {
        long long acc = 0;
        for (int j = c_NAC[i] - 1; j >= 0; j--) { long long c = c_AC[i][j] % l; if (c < 0) c += l; acc = (acc * tau + c) % l; }
        a[i] = acc;
    }
    const long long b2 = (a[0] * a[0] + 4 * a[1]) % l, b4 = (a[0] * a[2] + 2 * a[3]) % l, b6 = (a[2] * a[2] + 4 * a[4]) % l;     // a = [a1,a2,a3,a4,a6]
    const int c3 = 4 % l, c2 = (int)b2, c1 = (int)((2 * b4) % l), c0 = (int)b6;
    int f = c0, d1 = (c3 + c2 + c1) % l, d2 = (6 * c3 + 2 * c2) % l;
    const int d3 = (6 * c3) % l;
    int s = 0;
    for (int x = 0; x < l; x++) {
        if (f != 0) s += ((qr[f >> 5] >> (f & 31)) & 1) ? 1 : -1;
        f += d1;  if (f >= l)  f -= l;
        d1 += d2; if (d1 >= l) d1 -= l;
        d2 += d3; if (d2 >= l) d2 -= l;
    }
    out[tau] = (int16_t)(-s);                                // a_l = -(character sum)
}

typedef struct { int16_t p, q; float s; } ent_cmp_t;
int cmp_desc(const void *x, const void *y) { float a = ((const ent_cmp_t *)x)->s, b = ((const ent_cmp_t *)y)->s; return a < b ? 1 : (a > b ? -1 : 0); }
static int is_prime(int n) { if (n < 2) return 0; for (int d = 2; (long long)d * d <= n; d++) if (n % d == 0) return 0; return 1; }

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: tabx_build build FAM tabdir [plo phi] | tabx_build check FAM tabdir listfile\n"); return 1; }
    int fam = atoi(argv[2]);
    load_family(fam); FAM = fam;
    setenv("KUB_TABDIR", argv[3], 1);
    if (!strcmp(argv[1], "filter")) {          // tabx_build filter FAM tabdir infile THR4 [nthreads]  : keep the candidate lines "p q S300 S2000 S10000 H ..." with S10000 + dS >= THR4, append S(10^5)
        if (!tabx_load(fam)) { fprintf(stderr, "no tabx file\n"); return 1; }
        double thr4 = atof(argv[5]); int nt = argc > 6 ? atoi(argv[6]) : 16;
        FILE *f = fopen(argv[4], "r"); if (!f) { perror("infile"); return 1; }
        size_t cap = 1 << 20, n = 0; char **lines = (char **)malloc(cap * sizeof(char *)); char buf[1024];
        while (fgets(buf, sizeof buf, f)) { if (n == cap) { cap *= 2; lines = (char **)realloc(lines, cap * sizeof(char *)); } lines[n++] = strdup(buf); }
        fclose(f);
        char **outl = (char **)calloc(n, sizeof(char *));
        #pragma omp parallel for num_threads(nt) schedule(dynamic, 256)
        for (long i = 0; i < (long)n; i++) {
            long long p, q; double a, b, c, h;
            if (sscanf(lines[i], "%lld %lld %lf %lf %lf %lf", &p, &q, &a, &b, &c, &h) != 6) continue;
            double sx = c + tabx_sum(p, q);
            if (sx < thr4) continue;
            char o[1100]; size_t L = strlen(lines[i]); while (L && (lines[i][L-1] == '\n' || lines[i][L-1] == ' ')) L--;
            snprintf(o, sizeof o, "%.*s %.2f\n", (int)L, lines[i], sx);
            outl[i] = strdup(o);
        }
        size_t kept = 0;
        for (size_t i = 0; i < n; i++) if (outl[i]) { fputs(outl[i], stdout); kept++; }
        fprintf(stderr, "filter: %zu of %zu lines kept (S(1e5) >= %.1f)\n", kept, n, thr4);
        return 0;
    }
    if (!strcmp(argv[1], "list")) {            // tabx_build list FAM tabdir B nthreads keepfrac symflag : all coprime (p,q), 1<=q<=B, 0<|p|<=B, sorted by S(10^5) descending, top keepfrac printed as "p q S1e5"
        int B = atoi(argv[4]); int nt = argc > 5 ? atoi(argv[5]) : 12; double keep = argc > 6 ? atof(argv[6]) : 0.4; int sym = argc > 7 ? atoi(argv[7]) : 0;   // sym=1: p <= q/2 (Z/12: t <-> 1-t)
        if (!tabx_load(fam)) { fprintf(stderr, "no tabx file\n"); return 1; }
        build_tables(nt);
        typedef struct { int16_t p, q; float s; } ent_t;
        size_t cap = (size_t)(1.25 * (double)B * (double)B) + 1000; ent_t *E = (ent_t *)malloc(cap * sizeof(ent_t)); size_t n = 0;
        #pragma omp parallel for num_threads(nt) schedule(dynamic, 4)
        for (int q = 1; q <= B; q++) {
            ent_t *loc = (ent_t *)malloc((size_t)(2 * B + 1) * sizeof(ent_t)); int nl = 0;
            for (int p = -B; p <= B; p++) {
                if (p == 0) continue;
                if (sym == 1 && p > q / 2) continue;
                if (gcd64(p, q) != 1) continue;
                double a, b, c; sums_scalar(p, q, &a, &b, &c);
                loc[nl].p = (int16_t)p; loc[nl].q = (int16_t)q; loc[nl].s = (float)(c + tabx_sum(p, q)); nl++;
            }
            #pragma omp critical
            { if (n + nl <= cap) { memcpy(E + n, loc, (size_t)nl * sizeof(ent_t)); n += nl; } }
            free(loc);
        }
        int cmp_desc(const void *x, const void *y);
        qsort(E, n, sizeof(ent_t), cmp_desc);
        size_t out = (size_t)((double)n * keep);
        for (size_t i = 0; i < out; i++) printf("%d %d %.2f\n", E[i].p, E[i].q, E[i].s);
        fprintf(stderr, "list: %zu pairs, printed the best %zu\n", n, out);
        return 0;
    }
    if (!strcmp(argv[1], "check")) {
        if (!tabx_load(fam)) { fprintf(stderr, "no tabx file\n"); return 1; }
        FILE *f = fopen(argv[4], "r"); char line[512];
        while (fgets(line, sizeof line, f)) { long long p, q; if (sscanf(line, "%lld %lld", &p, &q) == 2) printf("%lld %lld dS=%.4f\n", p, q, tabx_sum(p, q)); }
        return 0;
    }
    int plo = argc > 4 ? atoi(argv[4]) : 10000, phi = argc > 5 ? atoi(argv[5]) : 100000;
    CK(cudaMemcpyToSymbol(c_AC, AC, sizeof AC)); CK(cudaMemcpyToSymbol(c_NAC, NAC, sizeof NAC));
    char fn[512]; snprintf(fn, sizeof fn, "%s/tabx_%d.bin", argv[3], fam);
    FILE *fo = fopen(fn, "wb"); if (!fo) { perror("open"); return 1; }
    uint32_t hdr[5] = { TABX_MAGIC, (uint32_t)fam, (uint32_t)plo, (uint32_t)phi, 0 };
    fwrite(hdr, 4, 5, fo);
    int16_t *d_out; CK(cudaMalloc(&d_out, sizeof(int16_t) * (size_t)phi));
    int16_t *h_out = (int16_t *)malloc(sizeof(int16_t) * (size_t)phi);
    int np = 0; double t0 = now();
    for (int l = plo + 1; l < phi; l++) {
        if (l % 2 == 0 || !is_prime(l)) continue;
        int nb = (l + 255) / 256;
        build_prime<<<nb, 256, ((l + 31) / 32) * sizeof(uint32_t)>>>(l, d_out);
        CK(cudaGetLastError());
        CK(cudaMemcpy(h_out, d_out, sizeof(int16_t) * (size_t)l, cudaMemcpyDeviceToHost));
        uint32_t ll = (uint32_t)l; fwrite(&ll, 4, 1, fo); fwrite(h_out, sizeof(int16_t), (size_t)l, fo);
        np++;
        if (np % 500 == 0) fprintf(stderr, "prime %d  (%d done)  %.0fs\n", l, np, now() - t0);
    }
    fseek(fo, 4 * 4, SEEK_SET); uint32_t n32 = np; fwrite(&n32, 4, 1, fo);
    fclose(fo);
    fprintf(stderr, "wrote %s: %d primes in (%d,%d)  %.0fs\n", fn, np, plo, phi, now() - t0);
    return 0;
}
#endif
