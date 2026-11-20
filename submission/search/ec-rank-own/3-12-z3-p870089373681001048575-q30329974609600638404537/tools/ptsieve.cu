// ptsieve.cu -- GPU sieve for rational points on 2-covering quartics y^2 = G(m,n) (hunt-points brief).
// Own code (ratpoints-style idea: one periodic bit pattern per small modulus, AND on 64-bit words, early exit).
//
// Input file, one quartic per line (written by pts.gp:writequartics):
//   qid mmax nmax a4 a3 a2 a1 a0 nint lo1 hi1 [lo2 hi2 ...]      (mmax < 0: search 0 <= m <= |mmax| only)
// G(m,n) = a4 m^4 + a3 m^3 n + a2 m^2 n^2 + a1 m n^3 + a0 n^4 (decimal, any size); search 1 <= n <= nmax,
// |m| <= mmax; [lo_i, hi_i] are the real intervals of x = m/n where G(x,1) >= 0 (inf allowed; a margin is added here).
// Output (stdout or --out): "qid m n" for every coprime (m,n) whose G(m,n) is a square modulo every sieve modulus
// (necessary condition; the host checks exactly with PARI). n = 0 is left to the host.
//
// For each modulus M (prime powers: 64, 27, 25, 49, 121, 169, primes 17..251) the host builds
//   T_M[r][o] = 64-bit word, bit j set iff G(o + j, r) mod M is a square mod M (indices mod M),
// so the word of row n starting at m is T_M[n mod M][m mod M]. Moduli are applied in increasing order of their
// pass fraction for this quartic; the loop stops when the word is zero.
//
// Thermal: every run starts search/gpusieve/watchdog.py (own flag + heartbeat under work/points, --tag points) and
// waits between launches while the pause flag exists, the heartbeat is stale, or work/gpusieve/BENCH-Z6 exists.
//
// Build: TMPDIR=$R/work/tmp nvcc -O3 -arch=sm_121 -std=c++17 -o bin/ptsieve ptsieve.cu
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include <chrono>
#include <thread>
#include <unistd.h>
#include <signal.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>

#define CK(x) do { cudaError_t e_ = (x); if (e_ != cudaSuccess) { fprintf(stderr, "CUDA %s at %s:%d\n", cudaGetErrorString(e_), __FILE__, __LINE__); exit(3); } } while (0)
typedef unsigned long long u64;
typedef unsigned int u32;
static const char *ROOT = ".";

#define MAXMOD 64
#define MAXINT 4
struct ModP { u32 M, minv, off, base, rowoff; };   // minv = floor(2^32 / M) (+1 rounding handled in fmod32)
__constant__ ModP c_mod[MAXMOD];
__constant__ double c_lo[MAXINT], c_hi[MAXINT];
__constant__ int c_nint, c_nmod;

__device__ __forceinline__ u32 fmod32(u32 x, u32 M, u32 minv) {
    u32 q = __umulhi(x, minv);
    u32 r = x - q * M;
    while (r >= M) r -= M;
    return r;
}
__device__ __forceinline__ u64 gcd64(u64 a, u64 b) { while (b) { u64 t = a % b; a = b; b = t; } return a; }

struct Hit { int qid; long long m; long long n; };

// one thread = one 64-bit word (64 consecutive m) of one row n
__global__ void k_sieve(const u64 *__restrict__ T, long long m0, u32 nwords, long long mmax, long long n0, u32 nrows,
                        int qid, Hit *hits, u32 *nhits, u32 maxhits, u64 *stat) {
    u64 gid = (u64)blockIdx.x * blockDim.x + threadIdx.x;
    u64 tot = (u64)nwords * nrows;
    if (gid >= tot) return;
    u32 w = (u32)(gid % nwords);
    long long n = n0 + (long long)(gid / nwords);
    long long ma = m0 + 64LL * w, mb = ma + 63;
    if (mb > mmax) mb = mmax;
    // real intervals
    bool ok = (c_nint == 0);
    double dn = (double)n;
    for (int i = 0; i < c_nint && !ok; i++) {
        double lo = c_lo[i] * dn, hi = c_hi[i] * dn;
        double slo = fabs(lo) * 1e-12 + 2.0, shi = fabs(hi) * 1e-12 + 2.0;
        if ((double)mb >= lo - slo && (double)ma <= hi + shi) ok = true;
    }
    if (!ok) return;
    u64 bits = ~0ULL;
    if (mb - ma < 63) bits = (1ULL << (mb - ma + 1)) - 1;
    u32 x = 64u * w;   // offset from m0; base = m0 mod M
    u32 un = (u32)n;   // n < 2^32 (checked on host)
    int k = 0;
    for (; k < c_nmod; k++) {
        ModP P = c_mod[k];
        u32 o = fmod32(x, P.M, P.minv) + P.base; if (o >= P.M) o -= P.M;
        u32 r = fmod32(un, P.M, P.minv);
        bits &= __ldg(&T[P.off + r * P.M + o]);
        if (!bits) break;
    }
    if (!bits) return;
    while (bits) {
        int j = __ffsll((long long)bits) - 1; bits &= bits - 1;
        long long m = ma + j;
        u64 am = m < 0 ? (u64)(-m) : (u64)m;
        if (gcd64(am, (u64)n) != 1) continue;
        u32 i = atomicAdd(nhits, 1u);
        if (i < maxhits) { hits[i].qid = qid; hits[i].m = m; hits[i].n = n; }
    }
}


// v2: block = (row n, chunk of WPT*256 words). Shared memory holds row (n mod M) of the first KS tables; per-modulus
// row base and start offset computed once per block. Rows that are empty modulo some modulus are skipped.
#define WPT 8
#define SHW 4096
__global__ void __launch_bounds__(256) k_sieve2(const u64 *__restrict__ T, const unsigned char *__restrict__ rowok, long long m0, u32 nwords,
                        long long mmax, long long n0, int qid, Hit *hits, u32 *nhits, u32 maxhits, int ks) {
    __shared__ u64 sT[SHW];
    __shared__ u32 sRow[MAXMOD], sO0[MAXMOD], sSh[MAXMOD];
    __shared__ int sDead;
    __shared__ long long sLo[MAXINT], sHi[MAXINT];
    long long n = n0 + blockIdx.y;
    u32 w0 = blockIdx.x * (256u * WPT);
    int t = threadIdx.x;
    if (t == 0) sDead = 0;
    if (t < c_nint) {   // integer m-range of interval t on this row (one FP64 evaluation per block)
        double dn = (double)n, lo = c_lo[t] * dn, hi = c_hi[t] * dn, lim = 9.0e18;
        lo -= fabs(lo) * 1e-12 + 2.0; hi += fabs(hi) * 1e-12 + 2.0;
        lo = fmax(lo, -lim); hi = fmin(hi, lim);
        sLo[t] = (long long)floor(lo); sHi[t] = (long long)ceil(hi);
    }
    __syncthreads();
    if (t < c_nmod) {
        ModP P = c_mod[t];
        u32 r = (u32)(n % (long long)P.M);
        sRow[t] = P.off + r * P.M;
        u64 x = (u64)P.base + 64ULL * w0;
        sO0[t] = (u32)(x % P.M);
        if (!rowok[P.rowoff + r]) sDead = 1;
    }
    __syncthreads();
    if (sDead) return;
    // copy shared rows: modulus k < ks goes to sSh[k]
    if (t == 0) { u32 o = 0; for (int k = 0; k < ks; k++) { sSh[k] = o; o += c_mod[k].M; } }
    __syncthreads();
    for (int k = 0; k < ks; k++) { u32 M = c_mod[k].M; for (u32 i = t; i < M; i += 256) sT[sSh[k] + i] = T[sRow[k] + i]; }
    __syncthreads();
    for (int i = 0; i < WPT; i++) {
        u32 dw = t + 256u * i, w = w0 + dw;
        if (w >= nwords) break;
        long long ma = m0 + 64LL * w, mb = ma + 63;
        if (mb > mmax) mb = mmax;
        bool ok = (c_nint == 0);
        for (int j = 0; j < c_nint; j++) if (mb >= sLo[j] && ma <= sHi[j]) ok = true;
        if (!ok) continue;
        u64 bits = ~0ULL;
        if (mb - ma < 63) bits = (1ULL << (mb - ma + 1)) - 1;
        u32 x = 64u * dw;
        int k = 0;
        for (; k < ks; k++) {
            ModP P = c_mod[k];
            u32 o = fmod32(x + sO0[k], P.M, P.minv);
            bits &= sT[sSh[k] + o];
            if (!bits) break;
        }
        if (bits) for (; k < c_nmod; k++) {
            ModP P = c_mod[k];
            u32 o = fmod32(x + sO0[k], P.M, P.minv);
            bits &= __ldg(&T[sRow[k] + o]);
            if (!bits) break;
        }
        while (bits) {
            int j = __ffsll((long long)bits) - 1; bits &= bits - 1;
            long long m = ma + j;
            u64 am = m < 0 ? (u64)(-m) : (u64)m;
            if (gcd64(am, (u64)n) != 1) continue;
            u32 h = atomicAdd(nhits, 1u);
            if (h < maxhits) { hits[h].qid = qid; hits[h].m = m; hits[h].n = n; }
        }
    }
}


__global__ void __launch_bounds__(256) k_sieve3(const u64 *__restrict__ T, const unsigned char *__restrict__ rowok, long long m0, u32 nwords,
                        long long mmax, long long n0, int qid, Hit *hits, u32 *nhits, u32 maxhits, int ks, int wpt) {
    __shared__ u64 sT[SHW];
    __shared__ u32 sRow[MAXMOD], sO0[MAXMOD], sSh[MAXMOD];
    __shared__ int sDead;
    __shared__ long long sLo[MAXINT], sHi[MAXINT];
    long long n = n0 + blockIdx.y;
    u32 w0 = blockIdx.x * (256u * wpt);
    int t = threadIdx.x;
    if (t == 0) sDead = 0;
    if (t < c_nint) {   // integer m-range of interval t on this row (one FP64 evaluation per block)
        double dn = (double)n, lo = c_lo[t] * dn, hi = c_hi[t] * dn, lim = 9.0e18;
        lo -= fabs(lo) * 1e-12 + 2.0; hi += fabs(hi) * 1e-12 + 2.0;
        lo = fmax(lo, -lim); hi = fmin(hi, lim);
        sLo[t] = (long long)floor(lo); sHi[t] = (long long)ceil(hi);
    }
    __syncthreads();
    if (t < c_nmod) {
        ModP P = c_mod[t];
        u32 r = (u32)(n % (long long)P.M);
        sRow[t] = P.off + r * P.M;
        u64 x = (u64)P.base + 64ULL * w0;
        sO0[t] = (u32)(x % P.M);
        if (!rowok[P.rowoff + r]) sDead = 1;
    }
    __syncthreads();
    if (sDead) return;
    // copy shared rows: modulus k < ks goes to sSh[k]
    if (t == 0) { u32 o = 0; for (int k = 0; k < ks; k++) { sSh[k] = o; o += c_mod[k].M; } }
    __syncthreads();
    for (int k = 0; k < ks; k++) { u32 M = c_mod[k].M; for (u32 i = t; i < M; i += 256) sT[sSh[k] + i] = T[sRow[k] + i]; }
    __syncthreads();
    // 4 words in flight per thread (independent loads overlap); words dw = t + 256*(i + j), j < 4
    for (int i = 0; i < wpt; i += 4) {
        u64 bits[4]; long long ma[4]; u32 x[4]; bool any = false;
        #pragma unroll
        for (int j = 0; j < 4; j++) {
            u32 dw = t + 256u * (i + j), w = w0 + dw;
            bits[j] = 0; ma[j] = m0 + 64LL * w; x[j] = 64u * dw;
            if (w >= nwords) continue;
            long long mb = ma[j] + 63; if (mb > mmax) mb = mmax;
            bool ok = (c_nint == 0);
            for (int q = 0; q < c_nint; q++) if (mb >= sLo[q] && ma[j] <= sHi[q]) ok = true;
            if (!ok) continue;
            bits[j] = (mb - ma[j] < 63) ? ((1ULL << (mb - ma[j] + 1)) - 1) : ~0ULL;
            any = true;
        }
        if (!any) continue;
        int k = 0;
        for (; k < ks; k++) {
            ModP P = c_mod[k]; u32 so = sO0[k], sb = sSh[k];
            #pragma unroll
            for (int j = 0; j < 4; j++) bits[j] &= sT[sb + fmod32(x[j] + so, P.M, P.minv)];
            if (!(bits[0] | bits[1] | bits[2] | bits[3])) break;
        }
        if (k == ks) for (; k < c_nmod; k++) {
            ModP P = c_mod[k]; u32 so = sO0[k], rb = sRow[k];
            #pragma unroll
            for (int j = 0; j < 4; j++) if (bits[j]) bits[j] &= __ldg(&T[rb + fmod32(x[j] + so, P.M, P.minv)]);
            if (!(bits[0] | bits[1] | bits[2] | bits[3])) break;
        }
        #pragma unroll
        for (int j = 0; j < 4; j++) {
            u64 bb = bits[j];
            while (bb) {
                int jj = __ffsll((long long)bb) - 1; bb &= bb - 1;
                long long m = ma[j] + jj;
                u64 am = m < 0 ? (u64)(-m) : (u64)m;
                if (gcd64(am, (u64)n) != 1) continue;
                u32 h = atomicAdd(nhits, 1u);
                if (h < maxhits) { hits[h].qid = qid; hits[h].m = m; hits[h].n = n; }
            }
        }
    }
}

// ---------------- table construction on the GPU ----------------
// thread = (modulus k, row r, column m) over the concatenated M^2 blocks
__device__ __forceinline__ int find_mod(const u32 *m2off, int NM, u32 g) { int lo = 0, hi = NM - 1; while (lo < hi) { int mid = (lo + hi + 1) / 2; if (m2off[mid] <= g) lo = mid; else hi = mid - 1; } return lo; }
__global__ void k_eval(unsigned char *bit, const u32 *co, const u32 *mlist, const u32 *m2off, const u32 *moff, const u32 *sqoff,
                       const unsigned char *sq, int NM, size_t tot, unsigned long long *cnt, unsigned char *rowok) {
    size_t g = (size_t)blockIdx.x * blockDim.x + threadIdx.x; if (g >= tot) return;
    int k = find_mod(m2off, NM, (u32)g); u32 M = mlist[k]; u32 loc = (u32)g - m2off[k]; u32 r = loc / M, m = loc % M;
    u32 pM = 2; while (M % pM) pM++;
    const u32 *a = co + 5 * k;
    u32 m2 = m * m % M, m3 = m2 * m % M, m4 = m3 * m % M, r2 = r * r % M, r3 = r2 * r % M, r4 = r3 * r % M;
    u32 v = (a[0] * m4 + a[1] * m3 % M * r + a[2] * m2 % M * r2 + a[3] * m % M * r3 + a[4] * r4) % M;
    unsigned char b = sq[sqoff[k] + v];
    if (r % pM == 0 && m % pM == 0) b = 0;
    bit[g] = b;
    if (b) { atomicAdd(&cnt[k], 1ULL); rowok[moff[k] + r] = 1; }
}
__global__ void k_expand(const unsigned char *bit, u64 *T, const u32 *mlist, const u32 *m2off, int NM, size_t tot) {
    size_t g = (size_t)blockIdx.x * blockDim.x + threadIdx.x; if (g >= tot) return;
    int k = find_mod(m2off, NM, (u32)g); u32 M = mlist[k]; u32 loc = (u32)g - m2off[k]; u32 r = loc / M, o = loc % M;
    const unsigned char *br = bit + m2off[k] + (size_t)r * M;
    u64 w = 0; u32 idx = o;
    for (int j = 0; j < 64; j++) { if (br[idx]) w |= 1ULL << j; if (++idx == M) idx = 0; }
    T[g] = w;
}

// ---------------- watchdog (same pattern as search/selmer/selmer_gpu.cu) ----------------
static pid_t g_wd = 0;
static std::string g_flag, g_hb;
static void wd_start() {
    std::string wdir = std::string(ROOT) + "/work/points";
    mkdir(wdir.c_str(), 0775);
    g_flag = wdir + "/pause"; g_hb = wdir + "/heartbeat";
    unlink(g_flag.c_str());
    pid_t me = getpid();
    g_wd = fork();
    if (g_wd == 0) {
        std::string pid = std::to_string(me);
        std::string log = std::string(ROOT) + "/logs/thermal.tsv";
        std::string wd = std::string(ROOT) + "/search/gpusieve/watchdog.py";
        execlp("python3", "python3", wd.c_str(), "--pid", pid.c_str(), "--flag", g_flag.c_str(), "--heartbeat", g_hb.c_str(),
               "--log", log.c_str(), "--tag", "points", (char *)NULL);
        _exit(127);
    }
    for (int i = 0; i < 100; i++) { struct stat st; if (stat(g_hb.c_str(), &st) == 0 && time(NULL) - st.st_mtime < 5) return; std::this_thread::sleep_for(std::chrono::milliseconds(100)); }
    fprintf(stderr, "watchdog did not start (no heartbeat)\n"); exit(4);
}
static void wd_stop() { if (g_wd > 0) { kill(g_wd, SIGTERM); waitpid(g_wd, NULL, 0); g_wd = 0; } }
static double g_paused = 0;
static void wd_gate() {
    auto t0 = std::chrono::steady_clock::now();
    while (1) {
        struct stat st; bool flag = stat(g_flag.c_str(), &st) == 0;
        bool stale = !(stat(g_hb.c_str(), &st) == 0 && time(NULL) - st.st_mtime <= 10);
        bool bench = stat((std::string(ROOT) + "/work/gpusieve/BENCH-Z6").c_str(), &st) == 0;
        if (!flag && !stale && !bench) break;
        std::this_thread::sleep_for(std::chrono::milliseconds(500));
    }
    g_paused += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}
static void read_gpu(double &w, double &t) {
    FILE *f = popen("nvidia-smi --query-gpu=power.draw,temperature.gpu --format=csv,noheader,nounits", "r");
    w = -1; t = -1; if (f) { if (fscanf(f, "%lf, %lf", &w, &t) != 2) { w = -1; t = -1; } pclose(f); }
}
static void sigh(int) { wd_stop(); _exit(2); }

// ---------------- host ----------------
static u32 modstr(const std::string &s, u32 M) {   // signed decimal mod M
    bool neg = false; size_t i = 0; if (s[0] == '-') { neg = true; i = 1; } else if (s[0] == '+') i = 1;
    u64 r = 0; for (; i < s.size(); i++) r = (r * 10 + (s[i] - '0')) % M;
    return neg ? (u32)((M - r) % M) : (u32)r;
}
static std::vector<u32> modulus_list() {
    std::vector<u32> L = {64, 27, 25, 49, 121, 169};
    for (u32 p = 17; p < 256; p++) { bool pr = true; for (u32 d = 2; d * d <= p; d++) if (p % d == 0) { pr = false; break; } if (pr) L.push_back(p); }
    return L;
}
static double now() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }

int main(int argc, char **argv) {
    std::string in, out; double target_ms = 150; int nmodmax = MAXMOD; bool nowd = false; u32 maxhits = 1u << 22; int ksmax = MAXMOD; int kver = 3; long long segarg = 1LL << 31;
    for (int i = 1; i < argc; i++) {
        std::string a = argv[i];
        if (a == "--out") out = argv[++i];
        else if (a == "--ms") target_ms = atof(argv[++i]);
        else if (a == "--nmod") nmodmax = atoi(argv[++i]);
        else if (a == "--no-watchdog-i-know") nowd = true;
        else if (a == "--v1") kver = 1;
        else if (a == "--v2") kver = 2;
        else if (a == "--seg") segarg = atoll(argv[++i]);
        else if (a == "--ks") ksmax = atoi(argv[++i]);
        else if (a[0] != '-' && in.empty()) in = a;
        else { fprintf(stderr, "usage: ptsieve quartics.txt [--out hits.txt] [--ms 150] [--nmod K]\n"); return 1; }
    }
    FILE *fin = fopen(in.c_str(), "r"); if (!fin) { perror(in.c_str()); return 1; }
    FILE *fo = out.empty() ? stdout : fopen(out.c_str(), "w"); if (!fo) { perror(out.c_str()); return 1; }
    signal(SIGINT, sigh); signal(SIGTERM, sigh);
    CK(cudaSetDeviceFlags(cudaDeviceScheduleBlockingSync));
    if (!nowd) wd_start();
    std::vector<u32> MODS = modulus_list();
    // static per-run modulus data
    int NM0 = (int)MODS.size(); std::vector<u32> h_m2off(NM0), h_moff(NM0), h_sqoff(NM0); size_t totM2 = 0, totM = 0, totSq = 0;
    for (int k = 0; k < NM0; k++) { h_m2off[k] = (u32)totM2; h_moff[k] = (u32)totM; h_sqoff[k] = (u32)totSq; totM2 += (size_t)MODS[k] * MODS[k]; totM += MODS[k]; totSq += MODS[k]; }
    std::vector<unsigned char> h_sq(totSq, 0);
    for (int k = 0; k < NM0; k++) { u32 M = MODS[k]; for (u64 y = 0; y < M; y++) h_sq[h_sqoff[k] + (y * y) % M] = 1; }
    u32 *d_mlist, *d_m2off, *d_moff, *d_sqoff, *d_co; unsigned char *d_sq, *d_bit; unsigned long long *d_cnt;
    CK(cudaMalloc(&d_mlist, 4 * NM0)); CK(cudaMalloc(&d_m2off, 4 * NM0)); CK(cudaMalloc(&d_moff, 4 * NM0)); CK(cudaMalloc(&d_sqoff, 4 * NM0));
    CK(cudaMalloc(&d_co, 4 * 5 * NM0)); CK(cudaMalloc(&d_sq, totSq)); CK(cudaMalloc(&d_bit, totM2)); CK(cudaMalloc(&d_cnt, 8 * NM0));
    CK(cudaMemcpy(d_mlist, MODS.data(), 4 * NM0, cudaMemcpyHostToDevice)); CK(cudaMemcpy(d_m2off, h_m2off.data(), 4 * NM0, cudaMemcpyHostToDevice));
    CK(cudaMemcpy(d_moff, h_moff.data(), 4 * NM0, cudaMemcpyHostToDevice)); CK(cudaMemcpy(d_sqoff, h_sqoff.data(), 4 * NM0, cudaMemcpyHostToDevice));
    CK(cudaMemcpy(d_sq, h_sq.data(), totSq, cudaMemcpyHostToDevice));
    u64 *d_T = nullptr; unsigned char *d_RO = nullptr;
    CK(cudaMalloc(&d_T, 8 * totM2)); CK(cudaMalloc(&d_RO, totM));
    Hit *d_hits; u32 *d_nh; u64 *d_stat;
    CK(cudaMalloc(&d_hits, sizeof(Hit) * maxhits)); CK(cudaMalloc(&d_nh, 4)); CK(cudaMalloc(&d_stat, 64));
    std::vector<Hit> hhits(maxhits);
    char line[1 << 16];
    double tall0 = now(); long double pts_all = 0; double kt_all = 0; double wsum = 0, wmax = 0, tmax = 0; int wn = 0;
    while (fgets(line, sizeof line, fin)) {
        if (line[0] == '#' || strlen(line) < 5) continue;
        std::vector<std::string> tok; { char *s = strtok(line, " \t\n"); while (s) { tok.push_back(s); s = strtok(NULL, " \t\n"); } }
        int qid = atoi(tok[0].c_str()); long long mmax = atoll(tok[1].c_str()), nmax = atoll(tok[2].c_str());
        long long mmin = -mmax; if (mmax < 0) { mmax = -mmax; mmin = 0; }   // negative mmax: half range 0..|mmax| (even quartic)
        std::string co[5] = {tok[3], tok[4], tok[5], tok[6], tok[7]};   // a4 a3 a2 a1 a0
        int nint = atoi(tok[8].c_str()); double lo[MAXINT], hi[MAXINT];
        if (nint > MAXINT) { fprintf(stderr, "too many intervals\n"); return 1; }
        for (int i = 0; i < nint; i++) { lo[i] = strtod(tok[9 + 2 * i].c_str(), NULL); hi[i] = strtod(tok[10 + 2 * i].c_str(), NULL); }
        if (nmax >= (1LL << 31) || mmax >= (1LL << 50)) { fprintf(stderr, "box too large\n"); return 1; }
        // tables, built on the GPU: k_eval (bits + counts + row flags), k_expand (64-bit words)
        int NM = (int)MODS.size();
        std::vector<u32> hco(NM * 5);
        for (int k = 0; k < NM; k++) for (int j = 0; j < 5; j++) hco[k * 5 + j] = modstr(co[j], MODS[k]);
        CK(cudaMemcpy(d_co, hco.data(), 4 * hco.size(), cudaMemcpyHostToDevice));
        CK(cudaMemset(d_cnt, 0, 8 * NM)); CK(cudaMemset(d_RO, 0, totM));
        k_eval<<<(u32)((totM2 + 255) / 256), 256>>>(d_bit, d_co, d_mlist, d_m2off, d_moff, d_sqoff, d_sq, NM, totM2, d_cnt, d_RO);
        k_expand<<<(u32)((totM2 + 255) / 256), 256>>>(d_bit, d_T, d_mlist, d_m2off, NM, totM2);
        CK(cudaDeviceSynchronize());
        std::vector<unsigned long long> hcnt(NM); CK(cudaMemcpy(hcnt.data(), d_cnt, 8 * NM, cudaMemcpyDeviceToHost));
        std::vector<int> ord; for (int k = 0; k < NM; k++) if ((double)hcnt[k] / ((double)MODS[k] * MODS[k]) < 0.999) ord.push_back(k);
        std::sort(ord.begin(), ord.end(), [&](int x, int y) { return (double)hcnt[x] / ((double)MODS[x] * MODS[x]) < (double)hcnt[y] / ((double)MODS[y] * MODS[y]); });
        if ((int)ord.size() > nmodmax) ord.resize(nmodmax);
        if (getenv("PTS_VERBOSE")) { fprintf(stderr, "# q %d moduli:", qid); for (int k = 0; k < (int)ord.size() && k < 16; k++) fprintf(stderr, " %u:%.3f", MODS[ord[k]], (double)hcnt[ord[k]] / ((double)MODS[ord[k]] * MODS[ord[k]])); fprintf(stderr, "\n"); }
        ModP P[MAXMOD]; double fprod = 1; int ks = 0; u32 shw = 0; size_t tsz = totM2;
        for (size_t k = 0; k < ord.size(); k++) {
            int q = ord[k]; u32 M = MODS[q];
            P[k].M = M; P[k].minv = (u32)((1ULL << 32) / M); P[k].off = h_m2off[q]; P[k].base = 0; P[k].rowoff = h_moff[q];
            if (ks == (int)k && shw + M <= SHW) { ks++; shw += M; }
            fprod *= (double)hcnt[q] / ((double)M * M);
        }

        int nm = (int)ord.size();
        CK(cudaMemcpyToSymbol(c_mod, P, sizeof(ModP) * nm)); CK(cudaMemcpyToSymbol(c_nmod, &nm, 4));
        CK(cudaMemcpyToSymbol(c_nint, &nint, 4)); if (nint) { CK(cudaMemcpyToSymbol(c_lo, lo, 8 * nint)); CK(cudaMemcpyToSymbol(c_hi, hi, 8 * nint)); }
        // m-range in segments of at most 2^31 values (32-bit word offsets inside a segment)
        const long long SEG = segarg;
        double rate_w = 0; u64 per = 1ULL << 26; u32 nh_tot = 0; double kt = 0; double t0q = now();
        CK(cudaMemset(d_nh, 0, 4));
        for (long long s0 = mmin; s0 <= mmax; s0 += SEG) {
            long long send = std::min(s0 + SEG - 1, mmax);
            u32 nwords = (u32)((send - s0 + 1 + 63) / 64);
            for (int k = 0; k < nm; k++) P[k].base = (u32)(((s0 % (long long)P[k].M) + P[k].M) % P[k].M);
            CK(cudaMemcpyToSymbol(c_mod, P, sizeof(ModP) * nm));
            long long n = 1;
            while (n <= nmax) {
                u64 rows = std::max<u64>(1, per / nwords); if ((long long)rows > nmax - n + 1) rows = nmax - n + 1;
                if (rows > 65535) rows = 65535;
                if (!nowd) wd_gate();
                u64 thr = rows * (u64)nwords; u32 blocks = (u32)((thr + 255) / 256);
                double t1 = now();
                if (kver == 3) { int wpt = (int)std::min<u64>(128, ((nwords + 255) / 256 + 3) / 4 * 4); dim3 g((nwords + 256 * wpt - 1) / (256 * wpt), (u32)rows);
                    k_sieve3<<<g, 256>>>(d_T, d_RO, s0, nwords, send, n, qid, d_hits, d_nh, maxhits, std::min(ks, ksmax), wpt); }
                else if (kver == 1) k_sieve<<<blocks, 256>>>(d_T, s0, nwords, send, n, (u32)rows, qid, d_hits, d_nh, maxhits, d_stat);
                else { dim3 g((nwords + 256 * WPT - 1) / (256 * WPT), (u32)rows); k_sieve2<<<g, 256>>>(d_T, d_RO, s0, nwords, send, n, qid, d_hits, d_nh, maxhits, std::min(ks, ksmax)); }
                CK(cudaDeviceSynchronize());
                double dt = now() - t1; kt += dt;
                rate_w = thr / dt;
                per = (u64)std::min(std::max(rate_w * target_ms / 1000.0, 1e6), 4e9);
                n += rows;
            }
        }
        CK(cudaMemcpy(&nh_tot, d_nh, 4, cudaMemcpyDeviceToHost));
        u32 nget = std::min(nh_tot, maxhits);
        if (nget) CK(cudaMemcpy(hhits.data(), d_hits, sizeof(Hit) * nget, cudaMemcpyDeviceToHost));
        for (u32 i = 0; i < nget; i++) fprintf(fo, "%d %lld %lld\n", hhits[i].qid, hhits[i].m, hhits[i].n);
        fflush(fo);
        static double last_smi = -1e9; static double w = -1, t = -1;
        if (now() - last_smi > 5) { read_gpu(w, t); last_smi = now(); }
        if (w > 0) { wsum += w; wn++; wmax = std::max(wmax, w); } tmax = std::max(tmax, t);
        long double pts = (long double)nmax * (mmax - mmin + 1);
        pts_all += pts; kt_all += kt;
        fprintf(stderr, "# q %d mmax %lld nmax %lld pts %.3Lg moduli %d tabMB %.1f passfrac %.2e hits %u%s kernel %.3fs wall %.3fs -> %.3Lg pts/s; GPU %.1f W %.0f C\n",
                qid, mmax, nmax, pts, nm, tsz * 8 / 1e6, fprod, nh_tot, nh_tot > maxhits ? " (OVERFLOW)" : "", kt, now() - t0q, pts / kt, w, t);
    }
    fprintf(stderr, "# total pts %.4Lg kernel %.2fs wall %.2fs paused %.1fs -> %.4Lg pts/s (kernel); GPU avg %.1f W max %.1f W, max %.0f C\n",
            pts_all, kt_all, now() - tall0, g_paused, kt_all > 0 ? pts_all / kt_all : 0.0L, wn ? wsum / wn : -1.0, wmax, tmax);
    if (fo != stdout) fclose(fo);
    wd_stop();
    return 0;
}
