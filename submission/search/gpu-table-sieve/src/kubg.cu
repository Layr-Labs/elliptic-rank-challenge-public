// kubg.cu - GPU stage 1 of the Mestre-Nagao table sieve (own code), stage 2/3 + exact checks on the CPU (reuses kub.c)
// build:  nvcc -O3 -arch=sm_121 -Xcompiler -O3,-march=native -o kubg kubg.cu -lpthread      (arch: the Spark's GB10; check nvidia-smi)
// usage:  kubg sieve FAM Hmax qlo qhi thr1 thr2 thr3 [cpu_threads]       (same output format as kub)
#include <cuda_runtime.h>
#include <string>
#include <vector>
#include <algorithm>
#define KUB_NO_MAIN
#include "kub.c"
#define TABX_NO_MAIN
#include "tabx.cu"          // extended tables (primes 10^4..10^5) for the host-side S(10^5) cut, see tabx.cu
static double THR4 = -1e18; static int USE_X = 0;      // env KUB_THR4: also require S(10^4) + dS(10^4..10^5) >= THR4 (needs tabx_<FAM>.bin in KUB_TABDIR)

#define NP1 61
#define CHUNK (1LL << 26)          // candidates per work item
#define MAXITEMS 65536
#define ITEM_MARGIN 4096          // free item slots required before a q is started (one q adds up to a few hundred items in KUB_PROG mode)
#define HITCAP (1 << 24)          // device hits per batch (24 bytes each); thin families need KUB_BATCHLOG <= 32 to stay below it
static int NOCOP = 0;
// FAM 44 / FAM 5 outer-value restriction ("small-prime-rich" values; FAM 5: |p| < q fundamental domain, role B uses the symmetry t -> -1/t).  W = {7,11,...,47}.  KUB_QRICH=t: only outer values q with at least t primes of W dividing q are swept (all inner values
// p in the shell).  KUB_ROLE=A: q = a (the larger of the Legendre pair), p = b in [1,q-1]  (the default fundamental domain);  KUB_ROLE=B: q = b (the smaller one), p = a in (q, inf):
// the family is symmetric in (a,b), so the same tables apply.  Measured on the S-selected rank>=7 curves of the Legendre family: "h_a >= 3 or h_b >= 3" keeps 55 % of them
// in 3.4 % of the pairs, "h >= 4" 16 % in 0.34 %.  Roles A and B are run as two separate sieves.
static int QRICH = 0, QRICHMAX = 0; static int ROLEB = 0;      // KUB_QRICHMAX=m (optional): also require at most m of the primes to divide q (complement passes)
static const int WP[12] = {7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47};
// KUB_PROG=n (with KUB_QRICH=KUB_QRICHMAX=h, FAM 44 / FAM 5): for the outer values with exactly h primes of W the inner values are NOT swept densely but only the multiples of the products m of n further
// primes of W (union of progressions p = m*k, duplicates removed per batch): the class "outer h, inner >= n" of the pairs whose small bad primes are split between the two factors.
struct Sub { int mask; int64_t m; };
static std::vector<Sub> SUBS[4];
static int PROG = 0;
static inline int64_t floordiv64(int64_t a, int64_t b) { int64_t q = a / b; if ((a % b != 0) && ((a < 0) != (b < 0))) q--; return q; }
static int64_t QMOD = 1, QRES = 0;   // KUB_QMOD=m:r  -> only q = r mod m (systematic subsample of the domain)

struct Item { int64_t q; int64_t plo; int64_t len; int64_t stride; };   // candidates p = plo + stride*u, u in [0,len)
struct Hit { int64_t p; int64_t q; int32_t s; int32_t pad; };

__constant__ int c_primes[NP1];
__constant__ int c_step[NP1];
__constant__ uint32_t c_magic[NP1];
__constant__ int c_ofs[NP1 + 1];
#define NP2MAX 320
__device__ int d_pr2[NP2MAX];             // odd primes < 2000 (all stage-2 primes, including the 61 stage-1 primes)
__device__ int d_off2[NP2MAX + 1];        // offsets of tab/inv per prime
__device__ int8_t *d_tab2;                 // tab_l[tau], concatenated
__device__ uint16_t *d_inv2;               // inv_l[r], concatenated
__device__ unsigned long long d_cnt1;      // number of stage-1 hits (before stage 2)
struct P2 { int l; int off; uint32_t magic; int pad; };
__device__ P2 d_p2[NP2MAX];                // per stage-2 prime: l, table offset, floor(2^32 / l)
static int NP2 = 0;

// Stage 2 for one candidate (p, q): sum of the tab_l[t mod l] over the stage-2 primes beyond the NP1 stage-1 primes, spread over the 32 lanes of a warp;
// the returned value is the full warp sum (identical in all lanes).
// sP2[k] = {l, table offset, floor(2^32/l), inverse of (q mod l) mod l or 0 when l | q (then the prime contributes nothing)}, built once per item in shared memory.
// No 64-bit division and no dependent global load before the table lookup (the old loop did q % l, p % l (both 64-bit) and a global inverse lookup
// per prime and hit); all iterations are independent so the table loads overlap.
// Early exit: after the first two iterations (the 64 stage-2 primes below ~700) the partial sum S(700) = s1 + sum is compared with c_tA (scaled units);
// a candidate below it cannot realistically reach the S(2000) threshold (it would need a gain of > 10 from the primes 700..2000, mean ~8 sd ~2.6 for rank-9
// curves) and the function returns a large negative number.  The host recomputes every surviving hit exactly, so the output semantics are unchanged.
#define NIT2 ((NP2MAX - NP1 + 31) / 32)
__constant__ int c_tA;
__device__ __forceinline__ int stage2_term(uint32_t ap, bool neg, const P2 &pk) {
    const uint32_t lq = (uint32_t)pk.l;
    uint32_t a = __umulhi(ap, pk.magic);
    int r = (int)(ap - a * lq); if (r >= pk.l) r -= pk.l;
    const int pm = (neg && r) ? pk.l - r : r;
    const uint32_t x = (uint32_t)pm * (uint32_t)pk.pad;
    uint32_t b = __umulhi(x, pk.magic);
    int tau = (int)(x - b * lq); if (tau >= pk.l) tau -= pk.l;
    const int v = d_tab2[pk.off + tau];
    return pk.pad ? v : 0;
}
__device__ __forceinline__ int warp_sum(int v) {
    #pragma unroll
    for (int o = 16; o > 0; o >>= 1) v += __shfl_xor_sync(0xffffffffu, v, o);
    return v;
}
__device__ __forceinline__ int stage2_sum(int64_t p, int lane, int np2, const P2 *sP2, int s1) {
    if (p > -(1LL << 31) && p < (1LL << 31)) {
        const uint32_t ap = (uint32_t)(p < 0 ? -p : p);
        const bool neg = p < 0;
        int sa = 0;
        #pragma unroll
        for (int j = 0; j < 2; j++) {
            const int k = NP1 + lane + 32 * j;
            if (k < np2) sa += stage2_term(ap, neg, sP2[k]);
        }
        sa = warp_sum(sa);
        if (s1 + sa < c_tA) return -100000;
        int sb = 0;
        #pragma unroll
        for (int j = 2; j < NIT2; j++) {
            const int k = NP1 + lane + 32 * j;
            if (k < np2) sb += stage2_term(ap, neg, sP2[k]);
        }
        return sa + warp_sum(sb);
    }
    int s2 = 0;
    for (int k = NP1 + lane; k < np2; k += 32) {
        const P2 pk = sP2[k];
        if (pk.pad == 0) continue;
        int pm = (int)(p % pk.l); if (pm < 0) pm += pk.l;
        const int tau = (int)(((uint32_t)pm * (uint32_t)pk.pad) % (uint32_t)pk.l);
        s2 += d_tab2[pk.off + tau];
    }
    return warp_sum(s2);
}

// per-item prologue of the stage-2 records (one 64-bit modulo per prime and item, instead of one per prime and hit)
__device__ __forceinline__ void stage2_prologue(int64_t q, int np2, P2 *sP2) {
    for (int k = NP1 + threadIdx.x; k < np2; k += blockDim.x) {
        P2 g = d_p2[k];
        const int qm = (int)(q % g.l);
        g.pad = qm ? (int)d_inv2[g.off + qm] : 0;
        sP2[k] = g;
    }
}

#define CK(x) do { cudaError_t e = (x); if (e != cudaSuccess) { fprintf(stderr, "CUDA error %s at %s:%d\n", cudaGetErrorString(e), __FILE__, __LINE__); exit(2); } } while (0)

__global__ void __launch_bounds__(256) stage1(const Item *items, const int8_t *d_pat, const int *d_patofs, int thr, int thr2i, int np2, Hit *hits, int *nhits, int cap) {
    __shared__ int8_t srow[9216];
    __shared__ int srlo[NP1];
    __shared__ int sstep[NP1];
    __shared__ P2 sP2[NP2MAX];          // stage 2: per-prime record incl. the inverse of q mod l (0 when l | q)
    __shared__ int shu[512];            // stage-1 hits of the current outer iteration: candidate index u
    __shared__ int shs[512];            // and their stage-1 sums
    __shared__ int shn;
    Item it = items[blockIdx.x];
    const int stride = (int)it.stride;
    int64_t q = it.q;
    for (int i = 0; i < NP1; i++) {
        int l = c_primes[i];
        int r = (int)(q % l);
        const int8_t *src = d_pat + d_patofs[i] + (size_t)r * l;
        for (int j = threadIdx.x; j < l; j += blockDim.x) srow[c_ofs[i] + j] = src[j];
        if (threadIdx.x == 0) { int64_t m = it.plo % l; if (m < 0) m += l; srlo[i] = (int)m; sstep[i] = (int)((256LL * stride) % l); }
    }
    stage2_prologue(q, np2, sP2);
    __syncthreads();
    for (int64_t u0 = 0; u0 < it.len; u0 += 4096) {
        int S[16];
#pragma unroll
        for (int j = 0; j < 16; j++) S[j] = 0;
        uint32_t x = (uint32_t)((u0 + threadIdx.x) * stride);
        for (int i = 0; i < NP1; i++) {
            int l = c_primes[i];
            uint32_t qq = __umulhi(x, c_magic[i]);
            int r = (int)(x - qq * (uint32_t)l); if (r >= l) r -= l;
            int idx = r + srlo[i]; if (idx >= l) idx -= l;
            const int8_t *R = srow + c_ofs[i];
            int d = sstep[i];
#pragma unroll
            for (int j = 0; j < 16; j++) {
                S[j] += R[idx];
                idx += d; if (idx >= l) idx -= l;
            }
        }
        if (threadIdx.x == 0) shn = 0;
        __syncthreads();
#pragma unroll
        for (int j = 0; j < 16; j++) {
            int64_t u = u0 + threadIdx.x + 256 * j;
            if (u < it.len && S[j] >= thr) {
                int pos = atomicAdd(&shn, 1);
                if (pos < 512) { shu[pos] = (int)u; shs[pos] = S[j]; }
            }
        }
        __syncthreads();
        int nh = shn < 512 ? shn : 512;
        if (threadIdx.x == 0 && nh > 0) atomicAdd(&d_cnt1, (unsigned long long)nh);
        // warp-cooperative stage 2: S2000 = stage-1 sum + the primes between 300 and 2000
        const int lane = threadIdx.x & 31, wid = threadIdx.x >> 5;
        for (int h = wid; h < nh; h += 8) {
            int64_t p = it.plo + (int64_t)stride * shu[h];
            const int s2 = stage2_sum(p, lane, np2, sP2, shs[h]);
            if (lane == 0 && shs[h] + s2 >= thr2i) {
                int pos = atomicAdd(nhits, 1);
                if (pos < cap) { hits[pos].p = p; hits[pos].q = q; hits[pos].s = shs[h]; hits[pos].pad = shs[h] + s2; }
            }
        }
        __syncthreads();
    }
}


// ---- candidate volume per class (for enrichment analysis): parity of (p,q), (p mod 3, q mod 3), bins of t = p/q
static double vol_par[4], vol_m3[9], vol_tb[14];
static int64_t flr(int64_t x, int64_t m) { return x >= 0 ? x / m : -((-x + m - 1) / m); }
static int64_t cnt_mod(int64_t a, int64_t b, int m, int r) { return flr(b - r, m) - flr(a - 1 - r, m); }
static void add_volume(int64_t q, int64_t a, int64_t b) {
    for (int r = 0; r < 2; r++) vol_par[r * 2 + (int)(q & 1)] += (double)cnt_mod(a, b, 2, r);
    for (int r = 0; r < 3; r++) vol_m3[r * 3 + (int)(q % 3)] += (double)cnt_mod(a, b, 3, r);
    for (int k = 0; k < 14; k++) {       // t bins: (-inf,-1.5), [-1.5,-1.25) ... [1.25,1.5), [1.5,inf)
        double e0 = k == 0 ? -1e18 : -1.5 + 0.25 * (k - 1), e1 = k == 13 ? 1e18 : -1.5 + 0.25 * k;
        int64_t lo = k == 0 ? a : (int64_t)ceil(e0 * (double)q), hi = k == 13 ? b : (int64_t)ceil(e1 * (double)q) - 1;
        if (lo < a) lo = a; if (hi > b) hi = b;
        if (hi >= lo) vol_tb[k] += (double)(hi - lo + 1);
    }
}


// ---------------------------------------------------------------------------------------------
// packed stage-1 kernel: 4 consecutive candidates per 32-bit word, 16-bit lanes (biased by 128), rows resampled per item
// so that consecutive candidates are consecutive row entries (the item stride is absorbed in the resampling).
// shared memory: for every prime l a table W_l[idx] (idx = 0..l+2) of 4-byte windows of the biased resampled row.
#define SUMW 8273               // sum over the 61 stage-1 primes of l words
#define BIAS 128
__constant__ int c_woff[NP1 + 1];
__constant__ int c_step1024[NP1];
__global__ void __launch_bounds__(256) stage1p(const Item *items, const int8_t *d_pat, const int *d_patofs, int thr, int thr2i, int np2, Hit *hits, int *nhits, int cap) {
    __shared__ uint32_t sW[SUMW];
    __shared__ P2 sP2[NP2MAX];
    __shared__ int shu[512];
    __shared__ int shs[512];
    __shared__ int shn;
    Item it = items[blockIdx.x];
    const int stride = (int)it.stride;
    int64_t q = it.q;
    // prologue: resampled biased row  rr[m] = R[(srlo + stride*m) mod l] + BIAS,  windows W[idx] = rr[idx..idx+3] (indices mod l)
    for (int i = 0; i < NP1; i++) {
        const int l = c_primes[i];
        const int r = (int)(q % l);
        const int8_t *src = d_pat + d_patofs[i] + (size_t)r * l;
        int64_t m0 = it.plo % l; if (m0 < 0) m0 += l;
        const int ms = (int)((stride % l));
        for (int idx = threadIdx.x; idx < l; idx += blockDim.x) {
            uint32_t w = 0;
            #pragma unroll
            for (int b = 0; b < 4; b++) {
                int m = (idx + 256 * b) % l;                   // position in the resampled row: the b-th candidate of the word is 256*b further
                int pos = (int)((m0 + (int64_t)ms * m) % l);   // position in the original pattern row
                uint32_t v = (uint32_t)((int)src[pos] + BIAS);
                w |= v << (8 * b);
            }
            sW[c_woff[i] + idx] = w;
        }
    }
    stage2_prologue(q, np2, sP2);
    __syncthreads();
    for (int64_t u0 = 0; u0 < it.len; u0 += 4096) {
        uint32_t acc0[4], acc1[4];
        #pragma unroll
        for (int j = 0; j < 4; j++) { acc0[j] = 0; acc1[j] = 0; }
        const uint32_t g0 = (uint32_t)(u0 + threadIdx.x);              // word j covers candidates g0 + 1024*j + 256*b, b = 0..3
        for (int i = 0; i < NP1; i++) {
            const int l = c_primes[i];
            uint32_t qq = __umulhi(g0, c_magic[i]);
            int idx = (int)(g0 - qq * (uint32_t)l); if (idx >= l) idx -= l;
            const int st = c_step1024[i];
            const uint32_t *Wl = sW + c_woff[i];
            #pragma unroll
            for (int j = 0; j < 4; j++) {
                uint32_t w = Wl[idx];
                acc0[j] += __byte_perm(w, 0, 0x5140);                     // lanes: byte0, byte1 (zero-extended)
                acc1[j] += __byte_perm(w, 0, 0x5342);                     // lanes: byte2, byte3
                idx += st; if (idx >= l) idx -= l;
            }
        }
        if (threadIdx.x == 0) shn = 0;
        __syncthreads();
        const int base = NP1 * BIAS;
        #pragma unroll
        for (int j = 0; j < 4; j++) {
            #pragma unroll
            for (int b = 0; b < 4; b++) {
                uint32_t a = (b < 2) ? acc0[j] : acc1[j];
                int s = (int)((b & 1) ? (a >> 16) : (a & 0xFFFFu)) - base;
                int64_t u = (int64_t)u0 + (int64_t)threadIdx.x + 1024 * j + 256 * b;
                if (u < it.len && s >= thr) {
                    int pos = atomicAdd(&shn, 1);
                    if (pos < 512) { shu[pos] = (int)u; shs[pos] = s; }
                }
            }
        }
        __syncthreads();
        int nh = shn < 512 ? shn : 512;
        if (threadIdx.x == 0 && nh > 0) atomicAdd(&d_cnt1, (unsigned long long)nh);
        const int lane = threadIdx.x & 31, wid = threadIdx.x >> 5;
        for (int h = wid; h < nh; h += 8) {
            int64_t p = it.plo + (int64_t)stride * shu[h];
            const int s2 = stage2_sum(p, lane, np2, sP2, shs[h]);
            if (lane == 0 && shs[h] + s2 >= thr2i) {
                int pos = atomicAdd(nhits, 1);
                if (pos < cap) { hits[pos].p = p; hits[pos].q = q; hits[pos].s = shs[h]; hits[pos].pad = shs[h] + s2; }
            }
        }
        __syncthreads();
    }
}


// SWAR variant of the packed kernel (see the comment inside)
#define BIAS2 12                // stage-1 table entries are within +-12 (8 * 2 ln(l)/sqrt(l) <= 11.8 for every prime l >= 3): biased entries in [0, 24]
#define SWB 10                  // primes per byte-lane block: 10 * 24 = 240 < 256
__global__ void __launch_bounds__(256) stage1s(const Item *items, const int8_t *d_pat, const int *d_patofs, int thr, int thr2i, int np2, Hit *hits, int *nhits, int cap) {
    __shared__ uint32_t sW[SUMW];
    __shared__ P2 sP2[NP2MAX];
    __shared__ int shu[512];
    __shared__ int shs[512];
    __shared__ int shn;
    Item it = items[blockIdx.x];
    const int stride = (int)it.stride;
    int64_t q = it.q;
    // prologue: resampled biased row  rr[m] = R[(srlo + stride*m) mod l] + BIAS,  windows W[idx] = rr[idx..idx+3] (indices mod l)
    for (int i = 0; i < NP1; i++) {
        const int l = c_primes[i];
        const int r = (int)(q % l);
        const int8_t *src = d_pat + d_patofs[i] + (size_t)r * l;
        int64_t m0 = it.plo % l; if (m0 < 0) m0 += l;
        const int ms = (int)((stride % l));
        for (int idx = threadIdx.x; idx < l; idx += blockDim.x) {
            uint32_t w = 0;
            #pragma unroll
            for (int b = 0; b < 4; b++) {
                int m = (idx + 256 * b) % l;                   // position in the resampled row: the b-th candidate of the word is 256*b further
                int pos = (int)((m0 + (int64_t)ms * m) % l);   // position in the original pattern row
                uint32_t v = (uint32_t)((int)src[pos] + BIAS2);
                w |= v << (8 * b);
            }
            sW[c_woff[i] + idx] = w;
        }
    }
    stage2_prologue(q, np2, sP2);
    __syncthreads();
    for (int64_t u0 = 0; u0 < it.len; u0 += 4096) {
        uint32_t acc0[4], acc1[4];
        #pragma unroll
        for (int j = 0; j < 4; j++) { acc0[j] = 0; acc1[j] = 0; }
        const uint32_t g0 = (uint32_t)(u0 + threadIdx.x);              // word j covers candidates g0 + 1024*j + 256*b, b = 0..3
        // SWAR: the byte lanes of up to SWB biased table words (each lane <= 2*BIAS2) are added directly (no cross-lane carry, sum <= 240 < 256);
        // only one unpack to the 16-bit accumulators per block of SWB primes instead of one per prime
        for (int i0 = 0; i0 < NP1; i0 += SWB) {
            uint32_t bs0 = 0, bs1 = 0, bs2 = 0, bs3 = 0;
            const int i1 = (i0 + SWB < NP1) ? i0 + SWB : NP1;
            for (int i = i0; i < i1; i++) {
                const int l = c_primes[i];
                uint32_t qq = __umulhi(g0, c_magic[i]);
                int idx = (int)(g0 - qq * (uint32_t)l); if (idx >= l) idx -= l;
                const int st = c_step1024[i];
                const uint32_t *Wl = sW + c_woff[i];
                bs0 += Wl[idx]; idx += st; if (idx >= l) idx -= l;
                bs1 += Wl[idx]; idx += st; if (idx >= l) idx -= l;
                bs2 += Wl[idx]; idx += st; if (idx >= l) idx -= l;
                bs3 += Wl[idx];
            }
            acc0[0] += __byte_perm(bs0, 0, 0x5140); acc1[0] += __byte_perm(bs0, 0, 0x5342);
            acc0[1] += __byte_perm(bs1, 0, 0x5140); acc1[1] += __byte_perm(bs1, 0, 0x5342);
            acc0[2] += __byte_perm(bs2, 0, 0x5140); acc1[2] += __byte_perm(bs2, 0, 0x5342);
            acc0[3] += __byte_perm(bs3, 0, 0x5140); acc1[3] += __byte_perm(bs3, 0, 0x5342);
        }
        if (threadIdx.x == 0) shn = 0;
        __syncthreads();
        const int base = NP1 * BIAS2;
        #pragma unroll
        for (int j = 0; j < 4; j++) {
            #pragma unroll
            for (int b = 0; b < 4; b++) {
                uint32_t a = (b < 2) ? acc0[j] : acc1[j];
                int s = (int)((b & 1) ? (a >> 16) : (a & 0xFFFFu)) - base;
                int64_t u = (int64_t)u0 + (int64_t)threadIdx.x + 1024 * j + 256 * b;
                if (u < it.len && s >= thr) {
                    int pos = atomicAdd(&shn, 1);
                    if (pos < 512) { shu[pos] = (int)u; shs[pos] = s; }
                }
            }
        }
        __syncthreads();
        int nh = shn < 512 ? shn : 512;
        if (threadIdx.x == 0 && nh > 0) atomicAdd(&d_cnt1, (unsigned long long)nh);
        const int lane = threadIdx.x & 31, wid = threadIdx.x >> 5;
        // stage 2 for four hits per warp at a time (h0, h0+8, h0+16, h0+24): the table loads and the warp reductions of the four hits overlap
        for (int h0 = wid; h0 < nh; h0 += 32) {
            int64_t pp[4]; bool ok[4]; int sh[4]; bool small = true;
            #pragma unroll
            for (int c = 0; c < 4; c++) {
                const int h = h0 + 8 * c;
                ok[c] = h < nh;
                sh[c] = ok[c] ? shs[h] : 0;
                pp[c] = ok[c] ? it.plo + (int64_t)stride * shu[h] : 0;
                if (!(pp[c] > -(1LL << 31) && pp[c] < (1LL << 31))) small = false;
            }
            if (small) {
                uint32_t ap[4]; bool ng[4]; int sa[4];
                #pragma unroll
                for (int c = 0; c < 4; c++) { ap[c] = (uint32_t)(pp[c] < 0 ? -pp[c] : pp[c]); ng[c] = pp[c] < 0; sa[c] = 0; }
                #pragma unroll
                for (int j = 0; j < 2; j++) {
                    const int k = NP1 + lane + 32 * j;
                    P2 pk; pk.l = 0; pk.off = 0; pk.magic = 0; pk.pad = 0;
                    if (k < np2) pk = sP2[k];
                    #pragma unroll
                    for (int c = 0; c < 4; c++) sa[c] += stage2_term(ap[c], ng[c], pk);
                }
                #pragma unroll
                for (int c = 0; c < 4; c++) sa[c] = warp_sum(sa[c]);
                #pragma unroll
                for (int c = 0; c < 4; c++) {
                    if (ok[c] && sh[c] + sa[c] >= c_tA) {
                        int sb = 0;
                        #pragma unroll
                        for (int j = 2; j < NIT2; j++) {
                            const int k = NP1 + lane + 32 * j;
                            if (k < np2) sb += stage2_term(ap[c], ng[c], sP2[k]);
                        }
                        const int s2 = sa[c] + warp_sum(sb);
                        if (lane == 0 && sh[c] + s2 >= thr2i) {
                            int pos = atomicAdd(nhits, 1);
                            if (pos < cap) { hits[pos].p = pp[c]; hits[pos].q = q; hits[pos].s = sh[c]; hits[pos].pad = sh[c] + s2; }
                        }
                    }
                }
            } else {
                for (int c = 0; c < 4; c++) {
                    if (!ok[c]) continue;
                    const int s2 = stage2_sum(pp[c], lane, np2, sP2, sh[c]);
                    if (lane == 0 && sh[c] + s2 >= thr2i) {
                        int pos = atomicAdd(nhits, 1);
                        if (pos < cap) { hits[pos].p = pp[c]; hits[pos].q = q; hits[pos].s = sh[c]; hits[pos].pad = sh[c] + s2; }
                    }
                }
            }
        }
        __syncthreads();
    }
}

typedef struct { Hit *h; int lo, hi; char *out; size_t olen, ocap; uint64_t c1, c2, c3; } hostjob_t;

static void *host_thr(void *arg) {
    hostjob_t *J = (hostjob_t *)arg;
    for (int k = J->lo; k < J->hi; k++) {
        int64_t p = J->h[k].p, q = J->h[k].q;
        if (p == 0 || gcd64(p, q) != 1) continue;
        double h = hnaive(FAM, (double)p, (double)q);
        if (h >= Hmax || h < Hmin) continue;
        J->c1++;
        double s1, s2, s3;
        if (!sums_staged(p, q, thr2, thr3, &s1, &s2, &s3)) continue;
        if (s1 < thr1 || s2 < thr2) continue;
        J->c2++;
        if (s3 < thr3) continue;
        double sx = 0;
        if (USE_X) { sx = s3 + tabx_sum(p, q); if (sx < THR4) continue; }
        J->c3++;
        char line[160];
        int n = USE_X ? snprintf(line, sizeof line, "%lld %lld %.2f %.2f %.2f %.4f %.2f\n", (long long)p, (long long)q, s1, s2, s3, h, sx)
                      : snprintf(line, sizeof line, "%lld %lld %.2f %.2f %.2f %.4f\n", (long long)p, (long long)q, s1, s2, s3, h);
        if (J->olen + n + 1 > J->ocap) { J->ocap = J->ocap ? J->ocap * 2 : 4096; J->out = (char *)realloc(J->out, J->ocap); }
        memcpy(J->out + J->olen, line, n); J->olen += n;
    }
    return NULL;
}

int main(int argc, char **argv) {
    if (argc < 9 || strcmp(argv[1], "sieve")) { fprintf(stderr, "usage: kubg sieve FAM Hmax qlo qhi thr1 thr2 thr3 [cpu_threads]\n"); return 1; }
    FAM = atoi(argv[2]);
    load_family(FAM);
    EARLY = getenv("KUB_EARLY") ? 1 : 0;
    if (getenv("KUB_THR4")) { THR4 = atof(getenv("KUB_THR4")); if (!tabx_load(FAM)) { fprintf(stderr, "KUB_THR4 set but tabx_%d.bin not found in KUB_TABDIR\n", FAM); return 3; } USE_X = 1; fprintf(stderr, "extended cut: S(10^5) >= %.2f (%d primes in the tabx file)\n", THR4, tabx_n); }      // host stage of the GPU driver: exact thresholds only (no intermediate cuts), hits are rare
    if (getenv("KUB_NOSYM")) SYM = 0;
    if (getenv("KUB_NOCOP")) NOCOP = 1;
    if (getenv("KUB_QRICH")) {
        QRICH = atoi(getenv("KUB_QRICH")); if (getenv("KUB_QRICHMAX")) QRICHMAX = atoi(getenv("KUB_QRICHMAX"));
        if (getenv("KUB_PROG")) {
            PROG = atoi(getenv("KUB_PROG"));
            if ((FAM != 44 && FAM != 5 && FAM != 66 && FAM != 67) || PROG < 1 || PROG > 3) { fprintf(stderr, "KUB_PROG needs FAM 44, 5, 66 or 67 and 1 <= n <= 3\n"); return 3; }
            for (int mask = 1; mask < 4096; mask++) { int sz = __builtin_popcount(mask); if (sz > 3) continue; int64_t m = 1; for (int j = 0; j < 12; j++) if (mask >> j & 1) m *= WP[j]; SUBS[sz].push_back({mask, m}); }
            fprintf(stderr, "progression mode: inner values divisible by %d further primes of W\n", PROG);
        }
        const char *role = getenv("KUB_ROLE"); ROLEB = (role && (role[0] == 'B' || role[0] == 'b')) ? 1 : 0;
        if (FAM != 44 && FAM != 5 && FAM != 66 && FAM != 67 && FAM != 88 && FAM != 26 && FAM != 261 && FAM != 265) { fprintf(stderr, "KUB_QRICH needs FAM 44, 5, 66, 67, 88, 26, 261 or 265\n"); return 3; }      // 66 / 67 = Z/6 with the outer value p resp. p+q (see fams.txt)
        fprintf(stderr, "outer-value restriction: >= %d of the primes 7..47 divide q; role %c (q = %s)\n", QRICH, ROLEB ? 'B' : 'A', ROLEB ? "the smaller number, |p| > q" : "the larger number, |p| < q");
    } else if (getenv("KUB_ROLE")) { const char *role = getenv("KUB_ROLE"); ROLEB = (role[0] == 'B' || role[0] == 'b') ? 1 : 0; if (ROLEB && FAM != 44 && FAM != 5) { fprintf(stderr, "KUB_ROLE=B needs FAM 44 or FAM 5\n"); return 3; } }
    if (getenv("KUB_QMOD")) { long long m = 1, r = 0; sscanf(getenv("KUB_QMOD"), "%lld:%lld", &m, &r); QMOD = m; QRES = r; }
    if (getenv("KUB_HMIN")) Hmin = atof(getenv("KUB_HMIN"));
    init_grid();
    Hmax = atof(argv[3]); qlo = atoll(argv[4]); qhi = atoll(argv[5]);
    thr1 = atof(argv[6]); thr2 = atof(argv[7]); thr3 = atof(argv[8]);
    { int64_t qa = q_auto_cap(Hmax); if (qhi > qa) qhi = qa; fprintf(stderr, "q range %lld..%lld (automatic cap %lld)\n", (long long)qlo, (long long)qhi, (long long)qa); }
    int nthr = argc > 9 ? atoi(argv[9]) : 8;
    uint64_t BATCHMAX = 1ULL << (getenv("KUB_BATCHLOG") ? atoi(getenv("KUB_BATCHLOG")) : 36);     // candidates per GPU batch
    double t0 = now();
    build_tables(nthr);
    fprintf(stderr, "tables ready in %.1fs: %d odd primes (stage1 %d)\n", now() - t0, nprimes, n1);
    if (n1 != NP1) { fprintf(stderr, "expected %d stage-1 primes, got %d\n", NP1, n1); return 1; }
    // device tables
    int h_ofs[NP1 + 1], h_patofs[NP1], h_step[NP1];
    uint32_t h_magic[NP1];
    size_t patbytes = 0; int ofs = 0;
    for (int i = 0; i < NP1; i++) { h_patofs[i] = (int)patbytes; patbytes += (size_t)primes[i] * primes[i]; h_ofs[i] = ofs; ofs += primes[i]; h_step[i] = 256 % primes[i]; h_magic[i] = (uint32_t)((1ULL << 32) / primes[i]); }
    h_ofs[NP1] = ofs;
    if (ofs > 9216) { fprintf(stderr, "shared row buffer too small: %d\n", ofs); return 1; }
    int8_t *hpat = (int8_t *)malloc(patbytes);
    for (int i = 0; i < NP1; i++) memcpy(hpat + h_patofs[i], pat[i], (size_t)primes[i] * primes[i]);
    int8_t *d_pat; int *d_patofs; Item *d_items; Hit *d_hits; int *d_nhits;
    CK(cudaMalloc(&d_pat, patbytes)); CK(cudaMemcpy(d_pat, hpat, patbytes, cudaMemcpyHostToDevice));
    CK(cudaMalloc(&d_patofs, sizeof h_patofs)); CK(cudaMemcpy(d_patofs, h_patofs, sizeof h_patofs, cudaMemcpyHostToDevice));
    CK(cudaMalloc(&d_items, MAXITEMS * sizeof(Item)));
    CK(cudaMalloc(&d_hits, (size_t)HITCAP * sizeof(Hit)));
    CK(cudaMalloc(&d_nhits, sizeof(int)));
    CK(cudaMemcpyToSymbol(c_primes, primes, sizeof(int) * NP1));
    CK(cudaMemcpyToSymbol(c_step, h_step, sizeof h_step));
    CK(cudaMemcpyToSymbol(c_magic, h_magic, sizeof h_magic));
    CK(cudaMemcpyToSymbol(c_ofs, h_ofs, sizeof h_ofs));
    // stage-2 tables (primes < 2000): tab, inverse, offsets
    NP2 = n2;
    if (NP2 > NP2MAX) { fprintf(stderr, "too many stage-2 primes\n"); return 1; }
    {
        int h_pr2[NP2MAX], h_off2[NP2MAX + 1]; size_t tot2 = 0;
        for (int i = 0; i < NP2; i++) { h_pr2[i] = primes[i]; h_off2[i] = (int)tot2; tot2 += primes[i]; }
        h_off2[NP2] = (int)tot2;
        int8_t *h_tab2 = (int8_t *)malloc(tot2); uint16_t *h_inv2 = (uint16_t *)malloc(tot2 * sizeof(uint16_t));
        for (int i = 0; i < NP2; i++) { memcpy(h_tab2 + h_off2[i], tab[i], primes[i]); memcpy(h_inv2 + h_off2[i], inv[i], primes[i] * sizeof(uint16_t)); }
        int8_t *dt2; uint16_t *di2;
        CK(cudaMalloc(&dt2, tot2)); CK(cudaMemcpy(dt2, h_tab2, tot2, cudaMemcpyHostToDevice));
        CK(cudaMalloc(&di2, tot2 * sizeof(uint16_t))); CK(cudaMemcpy(di2, h_inv2, tot2 * sizeof(uint16_t), cudaMemcpyHostToDevice));
        CK(cudaMemcpyToSymbol(d_tab2, &dt2, sizeof(dt2)));
        CK(cudaMemcpyToSymbol(d_inv2, &di2, sizeof(di2)));
        CK(cudaMemcpyToSymbol(d_pr2, h_pr2, sizeof(int) * NP2));
        CK(cudaMemcpyToSymbol(d_off2, h_off2, sizeof(int) * (NP2 + 1)));
        {
            P2 h_p2[NP2MAX]; memset(h_p2, 0, sizeof h_p2);
            for (int i = 0; i < NP2; i++) { h_p2[i].l = primes[i]; h_p2[i].off = h_off2[i]; h_p2[i].magic = (uint32_t)((1ULL << 32) / primes[i]); }
            CK(cudaMemcpyToSymbol(d_p2, h_p2, sizeof h_p2));
        }
        free(h_tab2); free(h_inv2);
    }
    {
        int h_woff[NP1 + 1], h_s1024[NP1]; int acc = 0;
        for (int i = 0; i < NP1; i++) { h_woff[i] = acc; acc += primes[i]; h_s1024[i] = 1024 % primes[i]; }
        h_woff[NP1] = acc;
        if (acc != SUMW) { fprintf(stderr, "SUMW mismatch %d vs %d\n", acc, SUMW); return 1; }
        CK(cudaMemcpyToSymbol(c_woff, h_woff, sizeof h_woff));
        CK(cudaMemcpyToSymbol(c_step1024, h_s1024, sizeof h_s1024));
    }
    int use_packed = getenv("KUB_NOPACK") ? 0 : 1;
    int use_swar = (use_packed && !getenv("KUB_NOSWAR")) ? 1 : 0;
    if (use_swar) {       // the SWAR kernel needs every stage-1 table entry within +-BIAS2
        int mx = 0;
        for (int i = 0; i < NP1; i++) for (int t = 0; t < primes[i]; t++) { int v = tab[i][t]; if (v < 0) v = -v; if (v > mx) mx = v; }
        if (mx > BIAS2) { fprintf(stderr, "stage-1 table entry %d exceeds BIAS2=%d: falling back to the packed kernel\n", mx, BIAS2); use_swar = 0; }
        else fprintf(stderr, "stage-1 table entries within +-%d\n", mx);
    }
    fprintf(stderr, "stage-1 kernel: %s\n", use_swar ? "packed + SWAR byte-lane blocks" : use_packed ? "packed (4 candidates per word)" : "scalar");
    int thr2i = (int)ceil(thr2 * SCALE - 1e-9);
    { double ta = getenv("KUB_TA") ? atof(getenv("KUB_TA")) : 14.0;          // early exit on S(700) (see stage2_sum); KUB_TA=-9999 disables
      int tai = ta < -1000 ? -(1 << 28) : (int)ceil(ta * SCALE - 1e-9);
      CK(cudaMemcpyToSymbol(c_tA, &tai, sizeof tai)); fprintf(stderr, "stage-2 early exit: S(700) >= %s\n", ta < -1000 ? "off" : (std::to_string(ta)).c_str()); }
    Hit *hhits = (Hit *)malloc((size_t)HITCAP * sizeof(Hit));
    Item *hitems = (Item *)malloc(MAXITEMS * sizeof(Item));
    int thr1i = (int)ceil(thr1 * SCALE - 1e-9);
    uint64_t total_cand = 0, tot1 = 0, tot2 = 0, tot3 = 0, overflow = 0, gpu_s1 = 0;
    double tstart = now(), tlast = tstart;
    int64_t q = qlo;
    int nitems = 0; uint64_t batch_cand = 0;
    int64_t ilo[40], ihi[40];
    for (;;) {
        // fill a batch
        while (q <= qhi && nitems + ITEM_MARGIN < MAXITEMS && batch_cand < BATCHMAX) {
            if (QMOD > 1) { int64_t adj = (((QRES - q) % QMOD) + QMOD) % QMOD; q += adj; if (q > qhi) break; }
            if (QRICH) { int hq = 0; for (int j = 0; j < 12; j++) hq += (q % WP[j] == 0); if (hq < QRICH || (QRICHMAX && hq > QRICHMAX)) { q++; continue; } }
            int nint = p_shell(q, Hmin, Hmax, ilo, ihi);
            if (SYM && nint && FAM == 5 && ROLEB) {      // Z/5 role B (outer value = |p| of the original pair, inner p' with |p'| > q'): two intervals per input interval
                int64_t tlo[40], thi[40]; int m = 0;
                for (int iv = 0; iv < nint; iv++) {
                    int64_t lo = ilo[iv], hi = ihi[iv];
                    int64_t h1 = hi < -q - 1 ? hi : -q - 1;                  // p' <= -q-1
                    if (lo <= h1 && m < 40) { tlo[m] = lo; thi[m] = h1; m++; }
                    int64_t l2 = lo > q + 1 ? lo : q + 1;                    // p' >= q+1
                    if (l2 <= hi && m < 40) { tlo[m] = l2; thi[m] = hi; m++; }
                }
                for (int i = 0; i < m; i++) { ilo[i] = tlo[i]; ihi[i] = thi[i]; }
                nint = m;
            } else
            if (SYM && nint) {
                int m = 0;
                for (int iv = 0; iv < nint; iv++) {
                    int64_t lo = ilo[iv], hi = ihi[iv];
                    if (FAM == 5) { if (lo < -q) lo = -q; if (hi > q) hi = q; }
                    else if (FAM == 42) { if (lo < 1) lo = 1; }
                else if (FAM == 8) { if (hi > (q >> 1)) hi = q >> 1; }          // Z/8: t and 1-t give the same curve -> fundamental domain p <= q/2
                    else if (FAM == 44) { if (ROLEB) { if (lo < q + 1) lo = q + 1; } else { if (lo < 1) lo = 1; if (hi > q - 1) hi = q - 1; } }   // Legendre Z/2xZ/4: fundamental domain 1 <= p=b < q=a (role A) or q=b < p=a (role B)
                    if (lo <= hi) { ilo[m] = lo; ihi[m] = hi; m++; }
                }
                nint = m;
            }
            for (int iv = 0; iv < nint; iv++) {
                int64_t lo = ilo[iv], hi = ihi[iv];
                add_volume(q, lo, hi);
                // classes of p that can be coprime to q at the primes 2 and 3 (stride-s progressions)
                int ncls = 1; int64_t s = 1; int64_t res[2] = {0, 0};
                if (!NOCOP) {
                    if (q % 6 == 0) { s = 6; ncls = 2; res[0] = 1; res[1] = 5; }
                    else if (q % 2 == 0) { s = 2; ncls = 1; res[0] = 1; }
                    else if (q % 3 == 0) { s = 3; ncls = 2; res[0] = 1; res[1] = 2; }
                }
                if (PROG) {
                    int qmask = 0; for (int j = 0; j < 12; j++) if (q % WP[j] == 0) qmask |= 1 << j;
                    for (const Sub &sb : SUBS[PROG]) {
                        if (sb.mask & qmask) continue;
                        const int64_t m = sb.m;
                        // FAM 66 (Z/6, outer = p, inner = q): the inner discriminant factor is p+q, i.e. the progressions are p' = -q (mod m) with an offset; no 2/3 coprimality classes then (the host drops gcd > 1)
                        const int64_t r0 = (FAM == 66) ? ((((-q) % m) + m) % m) : 0;
                        const int64_t klo = -floordiv64(-(lo - r0), m), khi = floordiv64(hi - r0, m);      // ceil((lo-r0)/m) .. floor((hi-r0)/m)
                        if (klo > khi) continue;
                        const int ncl = (FAM == 66) ? 1 : ncls; const int64_t ss = (FAM == 66) ? 1 : s;
                        for (int c = 0; c < ncl; c++) {
                            int64_t first = (FAM == 66) ? klo : klo + (((res[c] - klo) % ss) + ss) % ss;
                            if (first > khi) continue;
                            int64_t cnt = (khi - first) / ss + 1;
                            for (int64_t u = 0; u < cnt; u += CHUNK) {
                                int64_t len = cnt - u; if (len > CHUNK) len = CHUNK;
                                if (nitems < MAXITEMS) { hitems[nitems].q = q; hitems[nitems].plo = r0 + m * (first + ss * u); hitems[nitems].len = len; hitems[nitems].stride = m * ss; nitems++; batch_cand += (uint64_t)len; }
                            }
                        }
                    }
                    continue;
                }
                for (int c = 0; c < ncls; c++) {
                    int64_t first = lo + (((res[c] - lo) % s) + s) % s;
                    if (first > hi) continue;
                    int64_t cnt = (hi - first) / s + 1;
                    for (int64_t u = 0; u < cnt; u += CHUNK) {
                        int64_t len = cnt - u; if (len > CHUNK) len = CHUNK;
                        if (nitems < MAXITEMS) { hitems[nitems].q = q; hitems[nitems].plo = first + s * u; hitems[nitems].len = len; hitems[nitems].stride = s; nitems++; batch_cand += (uint64_t)len; }
                    }
                }
            }
            q++;
        }
        if (nitems == 0) break;
        CK(cudaMemcpy(d_items, hitems, nitems * sizeof(Item), cudaMemcpyHostToDevice));
        CK(cudaMemset(d_nhits, 0, sizeof(int)));
        { unsigned long long z = 0; CK(cudaMemcpyToSymbol(d_cnt1, &z, sizeof z)); }
        if (use_swar) stage1s<<<nitems, 256>>>(d_items, d_pat, d_patofs, thr1i, thr2i, NP2, d_hits, d_nhits, HITCAP);
        else if (use_packed) stage1p<<<nitems, 256>>>(d_items, d_pat, d_patofs, thr1i, thr2i, NP2, d_hits, d_nhits, HITCAP);
        else stage1<<<nitems, 256>>>(d_items, d_pat, d_patofs, thr1i, thr2i, NP2, d_hits, d_nhits, HITCAP);
        CK(cudaGetLastError());
        CK(cudaDeviceSynchronize());
        { unsigned long long c1 = 0; CK(cudaMemcpyFromSymbol(&c1, d_cnt1, sizeof c1)); gpu_s1 += c1; }
        int nh = 0; CK(cudaMemcpy(&nh, d_nhits, sizeof(int), cudaMemcpyDeviceToHost));
        if (nh > HITCAP) { overflow += (uint64_t)(nh - HITCAP); nh = HITCAP; }
        CK(cudaMemcpy(hhits, d_hits, (size_t)nh * sizeof(Hit), cudaMemcpyDeviceToHost));
        total_cand += batch_cand;
        if (PROG && nh > 1) {      // the same (p,q) can come from several progressions: keep one
            std::sort(hhits, hhits + nh, [](const Hit &x, const Hit &y) { return x.q != y.q ? x.q < y.q : x.p < y.p; });
            nh = (int)(std::unique(hhits, hhits + nh, [](const Hit &x, const Hit &y) { return x.q == y.q && x.p == y.p; }) - hhits);
        }
        // host-side checks in parallel
        pthread_t th[64]; hostjob_t jobs[64]; int nt = nthr > 64 ? 64 : nthr;
        for (int t = 0; t < nt; t++) {
            memset(&jobs[t], 0, sizeof(hostjob_t));
            jobs[t].h = hhits; jobs[t].lo = (int)((int64_t)nh * t / nt); jobs[t].hi = (int)((int64_t)nh * (t + 1) / nt);
            pthread_create(&th[t], NULL, host_thr, &jobs[t]);
        }
        for (int t = 0; t < nt; t++) {
            pthread_join(th[t], NULL);
            if (jobs[t].olen) { fwrite(jobs[t].out, 1, jobs[t].olen, stdout); }
            tot1 += jobs[t].c1; tot2 += jobs[t].c2; tot3 += jobs[t].c3; free(jobs[t].out);
        }
        fflush(stdout);
        double tn = now();
        if (tn - tlast > 20) { fprintf(stderr, "q<%lld  cand %.3g  rate %.3g/s  hits %llu\n", (long long)q, (double)total_cand, total_cand / (tn - tstart), (unsigned long long)tot1); tlast = tn; }
        nitems = 0; batch_cand = 0;
        if (q > qhi) break;
    }
    fprintf(stderr, "CLASSVOL par(p%%2,q%%2): %.6g %.6g %.6g %.6g\n", vol_par[0], vol_par[1], vol_par[2], vol_par[3]);
    fprintf(stderr, "CLASSVOL m3(p%%3,q%%3): "); for (int k = 0; k < 9; k++) fprintf(stderr, "%.6g ", vol_m3[k]); fprintf(stderr, "\n");
    fprintf(stderr, "CLASSVOL tbin: "); for (int k = 0; k < 14; k++) fprintf(stderr, "%.6g ", vol_tb[k]); fprintf(stderr, "\n");
    fprintf(stderr, "sieve %.1fs  rate %.3g candidates/s\n", now() - tstart, (double)total_cand / (now() - tstart));
    if (overflow) fprintf(stderr, "!!! HIT-BUFFER OVERFLOW: %llu hits lost - rerun with a smaller KUB_BATCHLOG or tighter thresholds\n", (unsigned long long)overflow);
    fprintf(stderr, "done: candidates %llu  after S300 %llu  after S2000 %llu  output %llu  (hit-buffer overflow %llu)\n",
            (unsigned long long)total_cand, (unsigned long long)gpu_s1, (unsigned long long)tot1, (unsigned long long)tot3, (unsigned long long)overflow);
    return 0;
}
