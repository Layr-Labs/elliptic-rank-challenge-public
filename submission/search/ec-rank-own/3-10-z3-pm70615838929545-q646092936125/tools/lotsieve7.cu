// lotsieve7 (build-lottery): sieve6 + --family z6k|z5|leg (height of the integral model, singular fibres) and dedup over
// n-subsets of any size n (row mask popcount).  Otherwise identical to search/gpusieve/src/sieve6.cu (git 3fbf3e5).
// sieve6: progression-row Mestre-Nagao sieve (64-bit outer values), built for the Kubert Z/6 like-for-like timing test
// (FILTER-ACCEPTANCE sec. 6c) and for class/progression sweeps of our own families.  Derived from sieve3.cu (dense SWAR
// prefilter for l < 2000, exact stages 2 and 3); see that file for the kernel layout.
//
//   S_L = sum_{5 <= l < L} -a_l log(l)/l,  a_l from the table (index t = pk * qk^-1 mod l, index l: qk = 0 mod l).
//
// Candidates come in ROWS (qk 64-bit, M) and SPANS (row, i0, len): candidate (row, i) is pk = M i, qk.  M must be coprime to
// qk (the progression pk = 0 mod M), |i| < 2^31, qk < 2^40.  Then t = pk/qk = i * rho_l with rho_l = M qk^-1 mod l, so a row
// of consecutive i is an ordinary table row:
//   dense stage-1 rows: P_l[r] with r = qk M^-1 mod l (P_l[r][s] = W(s r^-1) = W(i M / qk)); r = 0 is the class qk = 0
//     mod l (shared with sieve3); l | M gives the extra constant row r = l (t = 0 for every i);
//   dense 300 <= l < 2000: temp row t(s) = w(s rho) built in shared memory per block (rho = -1 encodes l | M);
//   exact stages: per warp and prime, rho = (M mod l)(qk mod l)^-1 and its Shoup constant floor(rho 2^32 / l) (table SHT),
//     t = |i| rho - umulhi(|i|, rho') l (+ sign), class qk = 0 mod l -> index l.
// gcd(pk, qk) = gcd(i, qk) (M coprime to qk) is tested exactly in stage 2.
// Records: SieveOut with p = i and q = row index; the host turns them into (pk, qk), applies the exact height test of the
// region (HLO <= H < HHI of the integral model) and deduplicates progressions (pk kept in the row whose triple is the three
// smallest W-primes dividing pk; see z6/z6enum.cpp).
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>
#include <algorithm>
#include <array>
#include <ctime>
#include <csignal>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <sys/time.h>
#include <cuda_runtime.h>
#include "gpusieve.h"

#define CK(x) do { cudaError_t e_ = (x); if (e_ != cudaSuccess) { \
  fprintf(stderr, "CUDA %s at %s:%d\n", cudaGetErrorString(e_), __FILE__, __LINE__); exit(2);} } while (0)

static const char* ROOT = ".";
#define MAXP1 72
#define MAXP2 1280
#define GMAX 40000
#define NT2 1024
#define FULL 0xffffffffu

struct RowD { long long q; long long M; int tri; int pad; };   // M 64-bit (PREFIX.rows64); old PREFIX.rows (int32 M) still read
struct Seg6 { int row, i0, len, pad; };

__constant__ int c_l[MAXP1], c_soff[MAXP1];
__constant__ long long c_poff[MAXP1];
__constant__ int c_np;
__constant__ int c2_l[MAXP2], c2_off[MAXP2], c2_ioff[MAXP2];
__constant__ unsigned int c2_m[MAXP2];

// x mod l for 0 <= x < 2^31, l < 2^15: m = floor(2^32/l)+1, the quotient estimate is exact or one too big.
__device__ __forceinline__ unsigned int fmod32(unsigned int x, unsigned int m, unsigned int l) {
  int r = (int)(x - __umulhi(x, m) * l); return r < 0 ? r + l : r; }
// q mod l for 0 <= q < 2^40
__device__ __forceinline__ unsigned int qmod(long long q, unsigned int m, unsigned int l) {   // 0 <= q < 2^62
  const unsigned int hi = (unsigned int)(q >> 40), mid = (unsigned int)((q >> 20) & 0xFFFFF), lo = (unsigned int)(q & 0xFFFFF), r20 = fmod32(1u << 20, m, l);
  return fmod32(fmod32(fmod32(hi, m, l) * r20 + mid, m, l) * r20 + lo, m, l); }
__device__ __forceinline__ unsigned int gcd_u(unsigned int a, unsigned int b) {
  if (!a) return b; if (!b) return a;
  int sh = __ffs(a | b) - 1; a >>= __ffs(a) - 1;
  while (b) { b >>= __ffs(b) - 1; if (a > b) { unsigned int t = a; a = b; b = t; } b -= a; }
  return a << sh; }
__device__ __forceinline__ bool coprime_iq(int i, long long q) {
  const unsigned int a = (unsigned int)abs(i); if (a == 0) return q == 1;
  return gcd_u(a, (unsigned int)((unsigned long long)q % a)) == 1; }

struct __align__(16) PP { int l; unsigned int m; int off; int ioff; };
__constant__ int c_cube, c_nogcd, c_swap, c_wide, c_z4;   // c_z4 (z4b): row u, inner s, t = u D^-1 s^-2 mod l
__device__ const signed char* g_QR;          // z4b: kronecker(x, l) at INV offsets
__device__ const unsigned short* g_DV;       // z4b: D^-1 mod l stored at each prime's INV offset (0 when l | D)
__device__ const unsigned short* g_INVP;     // z4b: the INV table (inverses mod l at the INV offsets)   // c_wide (z3w): row value |b| = (M << 64) | (uint64) q, multiplier 1   // c_swap (z3b): row value = |b|, inner i = a, t = b a^-3 mod l   // z3: row residue = (a mod l)^3; no coprimality filter (the host checks weighted primitivity)
__device__ __forceinline__ unsigned int rowres(long long q, unsigned int m, unsigned int l) {
  unsigned int x = qmod(q, m, l); if (c_cube && x) { const unsigned int x2 = fmod32(x * x, m, l); x = fmod32(x2 * x, m, l); } return x; }
__device__ __forceinline__ unsigned int mod64u(unsigned long long x, unsigned int m, unsigned int l) {   // any 64-bit x
  const unsigned int c32 = fmod32(0xFFFFFFFFu, m, l) + 1;
  return fmod32(fmod32((unsigned int)(x >> 32), m, l) * c32 + fmod32((unsigned int)x, m, l), m, l); }
// |b| mod l for the swap families: z3b row value q; z3w (q, M) = (low, high) 64-bit halves of |b| < 2^126
__device__ __forceinline__ unsigned int bres_of(long long q, long long M, unsigned int m, unsigned int l) {
  if (!c_wide) return qmod(q, m, l);
  const unsigned int c32 = fmod32(0xFFFFFFFFu, m, l) + 1, c64 = fmod32(c32 * c32, m, l);
  return fmod32(mod64u((unsigned long long)M, m, l) * c64 + mod64u((unsigned long long)q, m, l), m, l); }
// rho = M qk^-1 mod l (0 if l | M); 0xFFFF if l | qk (class infinity)
__device__ __forceinline__ unsigned int rho_of(long long q, long long M, const PP& pp, const unsigned short* __restrict__ INV) {
  const unsigned int qm = rowres(q, pp.m, pp.l); if (!qm) return 0xFFFFu;
  const unsigned int mm = qmod(M, pp.m, pp.l); if (!mm) return 0;
  return fmod32(mm * (unsigned int)INV[pp.ioff + qm], pp.m, pp.l); }
// z4b weight for inner s = ai at a prime: w = rho = u D^-1 mod l, 0xFFFE when l | D (additive, a_l = 0); l | s: a_l = kronecker(-u D, l)
// = kronecker(-rho, l) times the weight of a_l = +1 stored at base[off + l]; else table index t = rho s^-2
__device__ __forceinline__ int z4_w(unsigned int ai, unsigned int w, const PP& pp, const unsigned short* __restrict__ INV, const short* base) {
  if (w == 0xFFFEu) return 0;
  if (w == 0xFFFCu) return base[pp.off + pp.l];      // l | M: s = 0 mod l for every candidate, kronecker(-uD, l) = +1
  if (w == 0xFFFBu) return -base[pp.off + pp.l];     //                                          = -1
  const unsigned int s2 = fmod32(ai, pp.m, pp.l);
  if (!s2) { const unsigned int x = w ? pp.l - w : 0; return (int)g_QR[pp.ioff + x] * base[pp.off + pp.l]; }
  const unsigned int iv = INV[pp.ioff + s2], iv2 = fmod32(iv * iv, pp.m, pp.l); return base[pp.off + fmod32(w * iv2, pp.m, pp.l)]; }
// z4b row constant: rho = u D^-1 M^-2 mod l (inner s = M i); 0xFFFE: l | D; 0xFFFC / 0xFFFB: l | M with kronecker(-uD, l) = +1 / -1
__device__ __forceinline__ unsigned int z4_rho(long long u, long long M, int ioff, unsigned int m, unsigned int l) {
  const unsigned int d = g_DV[ioff]; if (!d) return 0xFFFEu;
  const unsigned int r0 = fmod32(qmod(u, m, l) * d, m, l), mm = qmod(M, m, l);
  if (!mm) { const unsigned int x = r0 ? l - r0 : 0; const int k = g_QR[ioff + x]; return k > 0 ? 0xFFFCu : 0xFFFBu; }
  const unsigned int mi = (unsigned int)g_INVP[ioff + mm], mi2 = fmod32(mi * mi, m, l); return fmod32(r0 * mi2, m, l); }
// z3b: index of t = bres * a^-3 mod l for inner a = +-ai (neg), class a = 0 mod l -> l.  z3w (c_wide): a = A0 + i with a0r = A0 mod l (true residue)
__device__ __forceinline__ unsigned int swap_t(unsigned int ai, bool neg, unsigned int bres, const PP& pp, const unsigned short* __restrict__ INV, unsigned int a0r) {
  unsigned int am = fmod32(ai, pp.m, pp.l);
  if (c_wide) { if (neg && am) am = pp.l - am; am += a0r; if (am >= (unsigned int)pp.l) am -= pp.l; if (!am) return pp.l;
    const unsigned int iv = INV[pp.ioff + am], iv3 = fmod32(fmod32(iv * iv, pp.m, pp.l) * iv, pp.m, pp.l); return fmod32(bres * iv3, pp.m, pp.l); }
  if (!am) return pp.l;
  const unsigned int iv = INV[pp.ioff + am], iv2 = fmod32(iv * iv, pp.m, pp.l), iv3 = fmod32(iv2 * iv, pp.m, pp.l);
  unsigned int t = fmod32(bres * iv3, pp.m, pp.l); return (neg && t) ? pp.l - t : t; }
// z3w row offset A0 (tri = low 32 bits, pad = high 32 bits) mod l as a true residue
__device__ __forceinline__ unsigned int a0res(const RowD& R, unsigned int m, unsigned int l) {
  const long long A0 = ((long long)R.pad << 32) | (long long)(unsigned int)R.tri; if (!A0) return 0;
  if (A0 > 0) return mod64u((unsigned long long)A0, m, l); const unsigned int r = mod64u((unsigned long long)(-A0), m, l); return r ? l - r : 0; }

// Sum of weights over the primes sp[0..n) for candidate (row, i).  uni: all lanes of the warp share row0.
__device__ __forceinline__ int groupsum6(const PP* sp, int n, unsigned int ai, bool neg, int row, int row0, bool uni,
    const short* base, const unsigned short* __restrict__ INV, const unsigned int* __restrict__ SHT, const RowD* __restrict__ rows, int lane) {
  int s = 0;
  if (uni) {
    const RowD R = rows[row0];
    for (int b = 0; b < n; b += 32) {
      unsigned int wj = 0xFFFFu, wpj = 0; const int ii = b + lane;
      if (ii < n) { const PP pp = sp[ii]; if (c_z4) wj = z4_rho(R.q, R.M, pp.ioff, pp.m, pp.l); else if (c_swap) { wj = bres_of(R.q, R.M, pp.m, pp.l); wpj = c_wide ? a0res(R, pp.m, pp.l) : 0; } else { wj = rho_of(R.q, R.M, pp, INV); if (wj != 0xFFFFu) wpj = SHT[pp.ioff + wj]; } }
      const int cnt = min(32, n - b);
      for (int u = 0; u < cnt; u++) {
        const unsigned int w = __shfl_sync(FULL, wj, u), wp = __shfl_sync(FULL, wpj, u);
        const PP pp = sp[b + u];
        unsigned int t;
        if (c_z4) { s += z4_w(ai, w, pp, INV, base); continue; }
        if (c_swap) t = swap_t(ai, neg, w, pp, INV, wp);
        else if (w != 0xFFFFu) { t = ai * w - __umulhi(ai, wp) * (unsigned int)pp.l; t = t >= (unsigned int)pp.l ? t - pp.l : t; t = (neg && t) ? pp.l - t : t; }
        else t = pp.l;
        s += base[pp.off + t];
      }
    }
  } else {
    const RowD R = rows[row];
    for (int k = 0; k < n; k++) {
      const PP pp = sp[k]; unsigned int t;
      if (c_z4) { s += z4_w(ai, z4_rho(R.q, R.M, pp.ioff, pp.m, pp.l), pp, INV, base); continue; }
      if (c_swap) { t = swap_t(ai, neg, bres_of(R.q, R.M, pp.m, pp.l), pp, INV, c_wide ? a0res(R, pp.m, pp.l) : 0); s += base[pp.off + t]; continue; }
      const unsigned int w = rho_of(R.q, R.M, pp, INV);
      if (w != 0xFFFFu) { const unsigned int wp = SHT[pp.ioff + w]; t = ai * w - __umulhi(ai, wp) * (unsigned int)pp.l; t = t >= (unsigned int)pp.l ? t - pp.l : t; t = (neg && t) ? pp.l - t : t; }
      else t = pp.l;
      s += base[pp.off + t];
    }
  }
  return s;
}

#define MAXGP 256
// Stage-2 group kernel.  mode 0: input (i,row), from scratch, gcd filter, cut, compact to (op,oq,os = S300).
// mode 1: iw[k] = is[k] + group sum.  mode 2: is[k] + group sum, cut, compact to SieveOut (S300 from os[k]).
__global__ void __launch_bounds__(NT2) stage2g(int mode, const int* __restrict__ ip, const int* __restrict__ iq, const int* __restrict__ is, int* __restrict__ iw,
    unsigned int n, int i0, int i1, int cut, const short* __restrict__ W2, const unsigned short* __restrict__ INV, const unsigned int* __restrict__ SHT,
    const RowD* __restrict__ rows, int* __restrict__ op, int* __restrict__ oq, int* __restrict__ os, SieveOut* __restrict__ orec,
    unsigned int* __restrict__ nout, unsigned int cap) {
  extern __shared__ short tab[];
  __shared__ PP sp[MAXGP];
  const int gbase = c2_off[i0], gsize = c2_off[i1 - 1] + c2_l[i1 - 1] + 1 - gbase, np = i1 - i0;
  for (int e = threadIdx.x; e < gsize; e += NT2) tab[e] = W2[gbase + e];
  for (int e = threadIdx.x; e < np; e += NT2) { PP pp; pp.l = c2_l[i0 + e]; pp.m = c2_m[i0 + e]; pp.off = c2_off[i0 + e] - gbase; pp.ioff = c2_ioff[i0 + e]; sp[e] = pp; }
  __syncthreads();
  const int lane = threadIdx.x & 31;
  const unsigned int nround = (n + gridDim.x * NT2 - 1) / (gridDim.x * NT2);
  for (unsigned int rd = 0; rd < nround; rd++) {
    const unsigned int k = (rd * gridDim.x + blockIdx.x) * NT2 + threadIdx.x;
    const bool valid = k < n;
    int i = 0, row = 0, s = 0;
    if (valid) { i = ip[k]; row = iq[k]; if (mode != 0) s = is[k]; }
    const int row0 = __shfl_sync(FULL, row, 0);
    const bool uni = __all_sync(FULL, !valid || row == row0);
    s += groupsum6(sp, np, (unsigned int)abs(i), i < 0, valid ? row : row0, row0, uni, tab, INV, SHT, rows, lane);
    bool sv = false;
    if (valid) {
      if (mode == 1) iw[k] = s;
      else if (mode == 2 || c_nogcd || coprime_iq(i, rows[row].q)) sv = s >= cut;
    }
    if (mode != 1) {
      const unsigned int m = __ballot_sync(FULL, sv);
      if (m) {
        const int leader = __ffs(m) - 1; unsigned int b0 = 0;
        if (lane == leader) b0 = atomicAdd(nout, __popc(m));
        b0 = __shfl_sync(FULL, b0, leader);
        if (sv) { unsigned int pos = b0 + __popc(m & ((1u << lane) - 1));
          if (pos < cap) { if (mode == 0) { op[pos] = i; oq[pos] = row; os[pos] = s; }
                           else { SieveOut r; r.p = i; r.q = row; r.s300 = os[k]; r.s2000 = s; r.s10000 = 0; r.bound = -1; orec[pos] = r; } } }
      }
    }
  }
}

#define MAXP3 1024
__global__ void __launch_bounds__(256) stage3(const SieveOut* __restrict__ in, unsigned int n, int i0, int i1, int cut,
    const short* __restrict__ W2, const unsigned short* __restrict__ INV, const unsigned int* __restrict__ SHT, const RowD* __restrict__ rows,
    SieveOut* __restrict__ out, unsigned int* __restrict__ nout, unsigned int cap) {
  __shared__ PP sp[MAXP3];
  const int np = i1 - i0;
  for (int e = threadIdx.x; e < np; e += 256) { PP pp; pp.l = c2_l[i0 + e]; pp.m = c2_m[i0 + e]; pp.off = c2_off[i0 + e]; pp.ioff = c2_ioff[i0 + e]; sp[e] = pp; }
  __syncthreads();
  const unsigned int k = blockIdx.x * 256 + threadIdx.x;
  const int lane = threadIdx.x & 31;
  const bool valid = k < n;
  SieveOut r; r.p = 0; r.q = 0; r.s2000 = 0;
  if (valid) r = in[k];
  const int row0 = __shfl_sync(FULL, r.q, 0);
  const bool uni = __all_sync(FULL, !valid || r.q == row0);
  const int s = r.s2000 + groupsum6(sp, np, (unsigned int)abs(r.p), r.p < 0, valid ? r.q : row0, row0, uni, W2, INV, SHT, rows, lane);
  r.s10000 = s;
  const bool sv = valid && s >= cut;
  const unsigned int m = __ballot_sync(FULL, sv);
  if (m) {
    const int leader = __ffs(m) - 1; unsigned int b0 = 0;
    if (lane == leader) b0 = atomicAdd(nout, __popc(m));
    b0 = __shfl_sync(FULL, b0, leader);
    if (sv) { unsigned int pos = b0 + __popc(m & ((1u << lane) - 1)); if (pos < cap) out[pos] = r; }
  }
}

// ---- Dense prefilter (stage 1 + 1b), one block = one Seg6 of <= SEGD candidates (see sieve3.cu dense12).
#define NTD 512
#define ITD 16
#define NCD 2
#define DD (NTD * ITD)
#define CHD (DD * 4)
#define SEGD (CHD * NCD)
#define MAXPB 256
#define MAXGB 64
__constant__ int cb_l[MAXPB], cb_woff[MAXPB], cb_goff[MAXPB], cb_toff[MAXPB], cb_D[MAXPB], cb_step[MAXPB], cb_ioff[MAXPB];
__constant__ unsigned short cb_wp[MAXPB], cb_wm[MAXPB];   // z4b: dense weights (with bias) of a_l = +1 and -1
__constant__ unsigned int cb_m[MAXPB];
__constant__ int cb_gs[MAXGB + 1], cb_ng, cb_np;
__constant__ int cd_step[MAXP1];
__constant__ unsigned int cd_m[MAXP1];
__constant__ int cd_ioff[MAXP1];
__global__ void __launch_bounds__(NTD, 1) dense12(const Seg6* __restrict__ segs, const RowD* __restrict__ rows, const unsigned long long* __restrict__ P,
    const unsigned short* __restrict__ WB, const unsigned short* __restrict__ INV, unsigned int cut1v, unsigned int cut2v, int gwords,
    int* __restrict__ op, int* __restrict__ oq, unsigned int* __restrict__ nout, unsigned int cap) {
  extern __shared__ unsigned long long shw[];
  __shared__ int sinv[MAXPB];
  __shared__ int srow[MAXP1]; __shared__ int ssh[MAXP1];
  const Seg6 sg = segs[blockIdx.x];
  const RowD R = rows[sg.row];
  const int len = sg.len, p0 = sg.i0, tid = threadIdx.x, lane = tid & 31;
  if (tid < c_np) {   // stage-1 row index r = qk M^-1 mod l (l | M: constant row l)
    const unsigned int l = c_l[tid], m = cd_m[tid], mm = c_wide ? 1u : qmod(R.M, m, l);
    const unsigned int qm = c_wide ? bres_of(R.q, R.M, m, l) : rowres(R.q, m, l);   // class qk = 0 mod l (row 0) takes precedence over l | M (row l), as in the exact stages
    srow[tid] = !qm ? 0 : (mm ? (int)fmod32(qm * (unsigned int)INV[cd_ioff[tid] + mm], m, l) : (int)l);
    if (c_z4) { const unsigned int rr = z4_rho(R.q, R.M, cd_ioff[tid], m, l); srow[tid] = rr == 0xFFFEu ? (int)l : (rr == 0xFFFCu ? (int)l + 1 : (rr == 0xFFFBu ? (int)l + 2 : (int)rr)); }
    ssh[tid] = c_wide ? (int)a0res(R, m, l) : 0;
  }
  for (int i = tid; i < cb_np; i += NTD) {   // rho = M qk^-1 mod l; 0: l | qk (class infinity); -1: l | M (t = 0)
    const unsigned int l = cb_l[i], m = cb_m[i], qm = c_wide ? bres_of(R.q, R.M, m, l) : rowres(R.q, m, l), mm = c_wide ? 1u : qmod(R.M, m, l);
    sinv[i] = c_swap ? (int)qm : (!qm ? 0 : (!mm ? -1 : (int)fmod32(mm * (unsigned int)INV[cb_ioff[i] + qm], m, l)));
    if (c_z4) { const unsigned int rr = z4_rho(R.q, R.M, cb_ioff[i], m, l); sinv[i] = rr == 0xFFFEu ? -3 : (rr == 0xFFFCu ? -4 : (rr == 0xFFFBu ? -5 : (int)rr)); }
  }
  __syncthreads();
  for (int i = 0; i < c_np; i++) {
    const int l = c_l[i];
    const unsigned long long* src = P + c_poff[i] + (long long)srow[i] * l;
    for (int s = tid; s < l; s += NTD) shw[c_soff[i] + s] = src[s];
  }
  __syncthreads();
  unsigned long long acc[NCD][ITD];
#pragma unroll
  for (int c = 0; c < NCD; c++)
#pragma unroll
    for (int it = 0; it < ITD; it++) acc[c][it] = 0;
  for (int i = 0; i < c_np; i++) {
    const int l = c_l[i], step = cd_step[i];
    const unsigned long long* row = shw + c_soff[i];
#pragma unroll
    for (int c = 0; c < NCD; c++) {
      int idx = (p0 + c * CHD + tid) % l; if (idx < 0) idx += l; idx += ssh[i]; if (idx >= l) idx -= l;
#pragma unroll
      for (int it = 0; it < ITD; it++) { acc[c][it] += row[idx]; idx += step; idx = idx >= l ? idx - l : idx; }
    }
  }
  unsigned int mask[NCD * ITD * 4 / 32];
#pragma unroll
  for (int w = 0; w < NCD * ITD * 4 / 32; w++) mask[w] = 0;
#pragma unroll
  for (int c = 0; c < NCD; c++)
#pragma unroll
    for (int it = 0; it < ITD; it++)
#pragma unroll
      for (int j = 0; j < 4; j++) {
        const int b = (c * ITD + it) * 4 + j;
        if (((unsigned int)(acc[c][it] >> (16 * j)) & 0xffffu) >= cut1v) mask[b >> 5] |= 1u << (b & 31);
      }
  unsigned short* tmp = (unsigned short*)(shw + gwords);
  for (int g = 0; g < cb_ng; g++) {
    const int i0 = cb_gs[g], i1 = cb_gs[g + 1];
    __syncthreads();
    for (int i = i0; i < i1; i++) {
      const int l = cb_l[i], inv = sinv[i]; const unsigned int m = cb_m[i];
      unsigned short* tr = tmp + cb_toff[i];
      const unsigned short* wb = WB + cb_woff[i];
      for (int s = tid; s < l; s += NTD) {
        unsigned short v;
        if (c_z4) {   // z4b temp row: -3 = l | D (a_l = 0: WB[l] holds the bias of a = 0); s = 0: kronecker(-rho) times w(+1); else WB[rho s^-2]
          if (inv == -3) v = wb[l];
          else if (inv == -4) v = cb_wp[i]; else if (inv == -5) v = cb_wm[i];
          else if (!s) { const unsigned int x = inv ? l - inv : 0; const int k = g_QR[cb_ioff[i] + x]; v = k > 0 ? cb_wp[i] : (k < 0 ? cb_wm[i] : wb[l]); }
          else { const unsigned int iv = INV[cb_ioff[i] + s], iv2 = fmod32(iv * iv, m, l); v = wb[fmod32((unsigned int)inv * iv2, m, (unsigned int)l)]; } }
        else if (c_swap) { int s2 = s + (c_wide ? (int)a0res(R, m, l) : 0); if (s2 >= l) s2 -= l;
          if (!s2) v = wb[l]; else { const unsigned int iv = INV[cb_ioff[i] + s2], iv3 = fmod32(fmod32(iv * iv, m, l) * iv, m, l); v = wb[fmod32((unsigned int)inv * iv3, m, (unsigned int)l)]; } }
        else if (inv > 0) v = wb[fmod32((unsigned int)s * (unsigned int)inv, m, (unsigned int)l)];
        else if (inv == 0) v = (s || c_cube) ? wb[l] : 0;   // z3: class a = 0 mod l is index l whatever b (no kill value)
        else v = wb[0];
        tr[s] = v;
      }
    }
    __syncthreads();
    for (int i = i0; i < i1; i++) {
      const int l = cb_l[i], D = cb_D[i];
      const unsigned short* tr = tmp + cb_toff[i];
      unsigned long long* gw = shw + cb_goff[i];
      for (int s = tid; s < l; s += NTD) {
        int s1 = s + D; s1 = s1 >= l ? s1 - l : s1; int s2 = s1 + D; s2 = s2 >= l ? s2 - l : s2; int s3 = s2 + D; s3 = s3 >= l ? s3 - l : s3;
        gw[s] = (unsigned long long)tr[s] | ((unsigned long long)tr[s1] << 16) | ((unsigned long long)tr[s2] << 32) | ((unsigned long long)tr[s3] << 48);
      }
    }
    __syncthreads();
    for (int i = i0; i < i1; i++) {
      const int l = cb_l[i], step = cb_step[i];
      const unsigned long long* row = shw + cb_goff[i];
#pragma unroll
      for (int c = 0; c < NCD; c++) {
        int idx = (p0 + c * CHD + tid) % l; if (idx < 0) idx += l;
#pragma unroll
        for (int it = 0; it < ITD; it++) { acc[c][it] += row[idx]; idx += step; idx = idx >= l ? idx - l : idx; }
      }
    }
  }
#pragma unroll
  for (int c = 0; c < NCD; c++)
#pragma unroll
    for (int it = 0; it < ITD; it++)
#pragma unroll
      for (int j = 0; j < 4; j++) {
        const int b = (c * ITD + it) * 4 + j;
        const int k = c * CHD + j * DD + it * NTD + tid;
        const bool sv = k < len && ((mask[b >> 5] >> (b & 31)) & 1u) && (((unsigned int)(acc[c][it] >> (16 * j)) & 0xffffu) >= cut2v);
        const unsigned int m = __ballot_sync(FULL, sv);
        if (m) {
          unsigned int b0 = 0; const int leader = __ffs(m) - 1;
          if (lane == leader) b0 = atomicAdd(nout, __popc(m));
          b0 = __shfl_sync(FULL, b0, leader);
          if (sv) { unsigned int pos = b0 + __popc(m & ((1u << lane) - 1)); if (pos < cap) { op[pos] = p0 + k; oq[pos] = sg.row; } }
        }
      }
}

// ---------------------------------------------------------------- host
static double now() { struct timeval tv; gettimeofday(&tv, 0); return tv.tv_sec + tv.tv_usec * 1e-6; }
static std::string utc() { char b[32]; time_t t = time(0); strftime(b, sizeof b, "%Y-%m-%dT%H:%M:%SZ", gmtime(&t)); return b; }
static bool exists(const std::string& f) { struct stat st; return stat(f.c_str(), &st) == 0; }
static double mtime_age(const std::string& f) { struct stat st; if (stat(f.c_str(), &st)) return 1e9; return now() - (st.st_mtim.tv_sec + st.st_mtim.tv_nsec * 1e-9); }
static std::string slurp(const std::string& f) { FILE* fp = fopen(f.c_str(), "rb"); if (!fp) { fprintf(stderr, "cannot open %s\n", f.c_str()); exit(1); }
  std::string s; char b[65536]; size_t n; while ((n = fread(b, 1, sizeof b, fp)) > 0) s.append(b, n); fclose(fp); return s; }
static std::vector<long> json_ints(const std::string& js, const char* key) {
  std::vector<long> v; std::string k = std::string("\"") + key + "\""; size_t i = js.find(k);
  if (i == std::string::npos) { fprintf(stderr, "key %s missing in al.json\n", key); exit(1); }
  i = js.find('[', i); size_t e = js.find(']', i); const char* p = js.c_str() + i + 1;
  while (p < js.c_str() + e) { char* q; long x = strtol(p, &q, 10); if (q == p) { p++; continue; } v.push_back(x); p = q; }
  return v; }
static long modinv(long a, long m) { long g = m, x = 0, x1 = 1, a1 = ((a % m) + m) % m; while (a1) { long t = g / a1; long tmp = g - t * a1; g = a1; a1 = tmp; tmp = x - t * x1; x = x1; x1 = tmp; } return ((x % m) + m) % m; }
static volatile sig_atomic_t g_stop = 0;
static void on_sig(int) { g_stop = 1; }
static pid_t g_wd = 0;
static void kill_wd() { if (g_wd > 0) { kill(g_wd, SIGTERM); waitpid(g_wd, 0, 0); g_wd = 0; } }
static void last_thermal(double& gc, double& w, double& mx) {
  gc = w = mx = NAN;
  FILE* fp = fopen((std::string(ROOT) + "/logs/thermal.tsv").c_str(), "rb"); if (!fp) return;
  fseek(fp, 0, SEEK_END); long sz = ftell(fp); long st = sz > 4096 ? sz - 4096 : 0; fseek(fp, st, SEEK_SET);
  std::string s(sz - st, 0); size_t n = fread(&s[0], 1, s.size(), fp); fclose(fp); s.resize(n);
  size_t e = s.find_last_of('\n', s.size() >= 2 ? s.size() - 2 : 0); std::string ln = s.substr(e == std::string::npos ? 0 : e + 1);
  std::vector<std::string> f; size_t a = 0; while (true) { size_t b = ln.find('\t', a); f.push_back(ln.substr(a, b - a)); if (b == std::string::npos) break; a = b + 1; }
  if (f.size() >= 6) { gc = atof(f[2].c_str()); w = atof(f[3].c_str()); mx = atof(f[f.size() - 2].c_str()); } }
// integral-model heights of the build-lottery families (binary128 long double on aarch64):
//   z6k [qk-2pk, -pk qk, -pk qk (qk-pk), 0, 0]; z5 [q-p, -p q, -p q^2, 0, 0]; leg [0, a^2+b^2, 0, a^2 b^2, 0]
static int g_fam = 0;   // 0 z6k, 1 z5, 2 leg, 3 z3, 4 z3b (z3 with rows = |b|, inner = a), 5 z3w, 7 z4b (rows u, inner s, t = u/(D s^2))
static long long g_dconst = 1;   // z4b: the squarefree D of the region
static long double Hfam(long long pk, long long qk) {
  long double p = pk, q = qk, a1, a2, a3, a4 = 0, a6 = 0;
  if (g_fam == 0) { a1 = q - 2 * p; a2 = -p * q; a3 = -p * q * (q - p); }
  else if (g_fam == 1) { a1 = q - p; a2 = -p * q; a3 = -p * q * q; }
  else { a1 = 0; a2 = p * p + q * q; a3 = 0; a4 = p * p * q * q; }
  long double b2 = a1 * a1 + 4 * a2, b4 = 2 * a4 + a1 * a3, b6 = a3 * a3 + 4 * a6;
  long double c4 = b2 * b2 - 24 * b4, c6 = -b2 * b2 * b2 + 36 * b2 * b4 - 216 * b6;
  return std::max(3 * logl(fabsl(c4)), 2 * logl(fabsl(c6))) / 12; }
static bool g_conic = false; static long long CFi[6];   // --conic: (pk,qk) = (F1(u,v), F2(u,v)), u = i, v = row q
static __int128 gcd128(__int128 a, __int128 b) { if (a < 0) a = -a; if (b < 0) b = -b; while (b) { __int128 t = a % b; a = b; b = t; } return a; }
static std::string s128(__int128 x) { if (x == 0) return "0"; bool ng = x < 0; if (ng) x = -x; std::string r; while (x) { r += char('0' + (int)(x % 10)); x /= 10; } if (ng) r += '-'; return std::string(r.rbegin(), r.rend()); }
static long double Hfam128(__int128 P, __int128 Q) {   // integral-model height, scaled to stay in range
  long double p = (long double)P, q = (long double)Q, m = fabsl(p) > fabsl(q) ? fabsl(p) : fabsl(q); p /= m; q /= m;
  long double a1, a2, a3, a4 = 0;
  if (g_fam == 0) { a1 = q - 2 * p; a2 = -p * q; a3 = -p * q * (q - p); } else if (g_fam == 1) { a1 = q - p; a2 = -p * q; a3 = -p * q * q; } else { a1 = 0; a2 = p * p + q * q; a3 = 0; a4 = p * p * q * q; }
  long double b2 = a1 * a1 + 4 * a2, b4 = 2 * a4 + a1 * a3, b6 = a3 * a3, c4 = b2 * b2 - 24 * b4, c6 = -b2 * b2 * b2 + 36 * b2 * b4 - 216 * b6;
  return logl(m) + std::max(3 * logl(fabsl(c4)), 2 * logl(fabsl(c6))) / 12; }
static bool singular(long long pk, long long qk) {
  if (g_fam == 0) return pk == 0 || pk == qk || qk == -8 * pk;
  if (g_fam == 1) return pk == 0;
  return pk == 0 || pk == qk || pk == -qk; }
static const int Wp[12] = {7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47};

int main(int argc, char** argv) {
  std::string table, family = "z6k", outdir, runid, spansf, pairsf, tag = "bench", wdflag = std::string(ROOT) + "/work/gpusieve/pause",
              hb = std::string(ROOT) + "/work/gpusieve/heartbeat", stopf = std::string(ROOT) + "/work/gpusieve/STOP", code = "?";
  int lim1 = 300, lim2 = 2000, lim3 = 10000;
  double K2 = 16384, cut1 = -1e9, cut2 = -1e9, cut3 = -1e9, duty = 1.0, batch = 4e8, dauto = 0, hlo = -1e9, hhi = 1e9;
  bool resume = false, nowd = false, dedup = true; long maxbatches = -1;
  for (int i = 1; i < argc; i++) {
    std::string a = argv[i]; auto nx = [&]() { if (i + 1 >= argc) { fprintf(stderr, "missing value for %s\n", a.c_str()); exit(1); } return std::string(argv[++i]); };
    if (a == "--table") table = nx(); else if (a == "--family") family = nx(); else if (a == "--tag") tag = nx();
    else if (a == "--spans") spansf = nx(); else if (a == "--pairs") pairsf = nx();
    else if (a == "--cut1") cut1 = atof(nx().c_str()); else if (a == "--cut2") cut2 = atof(nx().c_str()); else if (a == "--cut3") cut3 = atof(nx().c_str());
    else if (a == "--hshell") { hlo = atof(nx().c_str()); hhi = atof(nx().c_str()); }
    else if (a == "--duty") duty = atof(nx().c_str()); else if (a == "--duty-auto") dauto = atof(nx().c_str()); else if (a == "--batch") batch = atof(nx().c_str());
    else if (a == "--outdir") outdir = nx(); else if (a == "--runid") runid = nx(); else if (a == "--resume") resume = true; else if (a == "--code") code = nx();
    else if (a == "--no-dedup") dedup = false; else if (a == "--max-batches") maxbatches = atol(nx().c_str());
    else if (a == "--stop") stopf = nx();
    else if (a == "--no-watchdog-i-know") nowd = true;
    else if (a == "--dconst") g_dconst = atoll(nx().c_str());
    else if (a == "--conic") { std::string cs = nx(); g_conic = true; for (int k = 0; k < 6; k++) { CFi[k] = atoll(cs.c_str()); size_t c = cs.find(','); cs = c == std::string::npos ? "" : cs.substr(c + 1); } }
    else { fprintf(stderr, "unknown arg %s\n", a.c_str()); return 1; }
  }
  if (table.empty() || outdir.empty() || (spansf.empty() && pairsf.empty())) {
    fprintf(stderr, "usage: sieve6 --table DIR --outdir DIR (--spans PREFIX | --pairs FILE) [--cut1 S] [--cut2 S] [--cut3 S] [--hshell HLO HHI]\n"
                    "  [--duty D | --duty-auto T_C] [--batch N] [--resume] [--tag T] [--code C] [--no-dedup] [--max-batches N] [--stop FILE (default work/gpusieve/STOP)]\n"
                    "  spans: PREFIX.rows (int64 qk, int32 M, int32 tri), PREFIX.spans (int32 row, i0, len) from z6/z6enum\n"
                    "  pairs: lines 'pk qk M' (M | pk, gcd(M,qk) = 1): exact sums, all cuts off, output valsums.txt 'pk qk s300 s2000 s10000' (K2 units)\n"); return 1; }
  if (family == "z6k") g_fam = 0; else if (family == "z5") g_fam = 1; else if (family == "leg") g_fam = 2; else if (family == "z3") g_fam = 3; else if (family == "z3b") g_fam = 4; else if (family == "z3w") g_fam = 5; else if (family == "z4b") g_fam = 7;
  else { fprintf(stderr, "unknown --family %s (z6k z5 leg)\n", family.c_str()); return 1; }
  if (runid.empty()) runid = outdir.substr(outdir.find_last_of('/') + 1);
  mkdir(outdir.c_str(), 0775);
  const bool valmode = !pairsf.empty();
  if (valmode) { cut1 = cut2 = cut3 = -1e9; }

  // ---- table
  std::string js = slurp(table + "/al.json"), bin = slurp(table + "/al.bin");
  std::vector<long> primes = json_ints(js, "primes"), offs = json_ints(js, "offsets");
  const int16_t* al = (const int16_t*)bin.data();
  std::vector<int> I1, I2;
  for (size_t i = 0; i < primes.size(); i++) { if (primes[i] >= 5 && primes[i] < lim1) I1.push_back(i); if (primes[i] >= 5 && primes[i] < lim3) I2.push_back(i); }
  const int np1 = I1.size(), np2 = I2.size();
  if (np1 > MAXP1 || np2 > MAXP2) { fprintf(stderr, "too many primes\n"); return 1; }
  if (primes.back() < lim3 - 100) { fprintf(stderr, "table too short\n"); return 1; }
  for (int k = 0; k < np1; k++) if (I1[k] != I2[k]) { fprintf(stderr, "prime order\n"); return 1; }
  std::vector<int> l2(np2), off2(np2), ioff2(np2); std::vector<unsigned int> m2(np2); std::vector<short> W2; std::vector<unsigned short> INV; std::vector<unsigned int> SHT;
  int lv1 = 0, lv2 = 0;
  for (int k = 0; k < np2; k++) {
    int i = I2[k]; int l = primes[i]; l2[k] = l; m2[k] = (unsigned int)((1ULL << 32) / l + 1); off2[k] = W2.size(); ioff2[k] = INV.size();
    for (int t = 0; t <= l; t++) { const double av = (g_fam == 7 && t == l) ? 1.0 : (double)al[offs[i] + t];   // z4b: slot l = weight of a_l = +1 (kronecker rule for l | s)
      long w = lround(K2 * (-av) * log((double)l) / l); if (labs(w) > 32767) { fprintf(stderr, "K2 too large\n"); return 1; } W2.push_back((short)w); }
    for (int r = 0; r < l; r++) { INV.push_back((unsigned short)(r ? modinv(r, l) : 0)); SHT.push_back((unsigned int)(((unsigned long long)r << 32) / l)); }
    if (l < lim1) lv1 = k + 1; if (l < lim2) lv2 = k + 1;
  }
  std::vector<std::pair<int, int>> groups; groups.push_back({0, lv1});
  { int a = lv1; while (a < lv2) { int b = a, sz = 0; while (b < lv2 && sz + l2[b] + 1 <= GMAX) { sz += l2[b] + 1; b++; } if (b == a) b++; groups.push_back({a, b}); a = b; } }
  for (auto& g : groups) if (g.second - g.first > MAXGP) { fprintf(stderr, "group too large\n"); return 1; }
  if (np2 - lv2 > MAXP3) { fprintf(stderr, "too many stage-3 primes\n"); return 1; }
  int maxg = 0; for (auto& g : groups) maxg = std::max(maxg, off2[g.second - 1] + l2[g.second - 1] + 1 - off2[g.first]);
  auto icut = [&](double c) { return c <= -1e8 ? INT32_MIN : (int)ceil(c * K2 - 1e-9); };
  const int cut1i = icut(cut1), cut2i = icut(cut2), cut3i = icut(cut3);

  // ---- dense prefilter tables (common scale KC for all primes < lim2); stage-1 rows: l+1 rows per prime (row l: t = 0)
  const int npb = lv2 - lv1;
  if (npb > MAXPB - 1) { fprintf(stderr, "dense: too many primes\n"); return 1; }
  auto mlv = [&](int i) { long l = primes[i]; int mx = 0; for (long t = 0; t <= l; t++) mx = std::max(mx, abs((int)al[offs[i] + t])); return mx * log((double)l) / l; };
  double sm = 0; for (int k = 0; k < np1; k++) sm += mlv(I1[k]); for (int k = lv1; k < lv2; k++) sm += mlv(I2[k]);
  const double KC = floor((65535.0 - np1 - npb) / (2 * sm));
  auto wc = [&](int i, int a) { long l = primes[i]; return (int)lround(KC * (-(double)a) * log((double)l) / l); };
  auto Bof = [&](int i) { long l = primes[i]; int mx = 0; for (long t = 0; t <= l; t++) mx = std::max(mx, abs(wc(i, al[offs[i] + t]))); return mx; };
  int sumBC1 = 0, sumBC = 0;
  std::vector<int> BC1(np1); for (int k = 0; k < np1; k++) { BC1[k] = Bof(I1[k]); sumBC1 += BC1[k]; }
  sumBC = sumBC1; std::vector<int> BCb(npb); for (int k = 0; k < npb; k++) { BCb[k] = Bof(I2[lv1 + k]); sumBC += BCb[k]; }
  if (2 * sumBC > 65535) { fprintf(stderr, "dense: lane overflow possible\n"); return 1; }
  const double errC1 = np1 * 0.5 / KC, errC2 = (np1 + npb) * 0.5 / KC;
  std::vector<int> L1(np1), soff(np1), dstep(np1), dioff(np1); std::vector<unsigned int> dm(np1); std::vector<long long> poff(np1); long long ptot = 0; int tot = 0;
  const int PCR = g_fam == 7 ? 3 : 1;   // z4b: rows l+1, l+2 = constant kronecker +1 / -1 (l | M)
  for (int k = 0; k < np1; k++) { L1[k] = primes[I1[k]]; soff[k] = tot; tot += L1[k]; poff[k] = ptot; ptot += (long long)L1[k] * (L1[k] + PCR);
    dstep[k] = NTD % L1[k]; dm[k] = m2[k]; dioff[k] = ioff2[k]; }
  std::vector<unsigned long long> PC(ptot);
  for (int k = 0; k < np1; k++) {
    int l = L1[k], D = (NTD * ITD) % l; std::vector<int> wb(l + 1), row(l);
    for (int t = 0; t <= l; t++) wb[t] = wc(I1[k], al[offs[I1[k]] + t]) + BC1[k];
    for (int r = 0; r <= l + PCR - 1; r++) { long ri = (r && r < l) ? modinv(r, l) : 0;
      if (g_fam == 7) { const int wpl = wc(I1[k], 1) + BC1[k], wmi = wc(I1[k], -1) + BC1[k];
        for (int s2 = 0; s2 < l; s2++) { if (r == l) { row[s2] = wb[l]; continue; }   // row l: l | D, a_l = 0 (wb[l] = bias: table class l is 0)
          if (r == l + 1) { row[s2] = wpl; continue; } if (r == l + 2) { row[s2] = wmi; continue; }   // l | M: constant kronecker weight
          if (s2) { long iv = modinv(s2, l); row[s2] = wb[(long)r * (iv * iv % l) % l]; }
          else if (r) { long x = (l - r) % l, kq = 1; for (long e = (l - 1) / 2, bb = x; e; e >>= 1, bb = bb * bb % l) if (e & 1) kq = kq * bb % l; row[s2] = kq == 1 ? wpl : wmi; }
          else row[s2] = 0; } }
      else if (g_fam == 4 || g_fam == 5) { for (int s2 = 0; s2 < l; s2++) { long iv = s2 ? modinv(s2, l) : 0; row[s2] = r == l ? wb[0] : (s2 ? wb[(long)r * (iv * iv % l * iv % l) % l] : wb[l]); } }
      else for (int s2 = 0; s2 < l; s2++) row[s2] = r == l ? wb[0] : (r ? wb[(s2 * ri) % l] : ((s2 || g_fam == 3) ? wb[l] : 0));
      for (int s2 = 0; s2 < l; s2++) { unsigned long long w = 0; for (int j = 0; j < 4; j++) w |= (unsigned long long)row[(s2 + j * D) % l] << (16 * j); PC[poff[k] + (long long)r * l + s2] = w; } }
  }
  std::vector<int> bl(npb), bwoff(npb), bgoff(npb), btoff(npb), bD(npb), bstep(npb), bioff(npb), bgs; std::vector<unsigned int> bm(npb); std::vector<unsigned short> WB;
  int maxsh = 0; CK(cudaDeviceGetAttribute(&maxsh, cudaDevAttrMaxSharedMemoryPerBlockOptin, 0));
  const int gw = std::min(10200, (maxsh - (int)sizeof(int) * (MAXPB + MAXP1) - 512) / 10);
  { int words = 0; bgs.push_back(0);
    for (int k = 0; k < npb; k++) {
      int l = primes[I2[lv1 + k]];
      if (words + l > gw && words > 0) { bgs.push_back(k); words = 0; }
      bl[k] = l; bgoff[k] = words; words += l; bD[k] = (NTD * ITD) % l; bstep[k] = NTD % l; bm[k] = m2[lv1 + k]; bioff[k] = ioff2[lv1 + k];
      bwoff[k] = WB.size(); for (int t = 0; t <= l; t++) WB.push_back((unsigned short)(wc(I2[lv1 + k], al[offs[I2[lv1 + k]] + t]) + BCb[k])); }
    bgs.push_back(npb); }
  if ((int)bgs.size() - 1 > MAXGB) { fprintf(stderr, "dense: too many groups\n"); return 1; }
  for (size_t g = 0; g + 1 < bgs.size(); g++) { int off = 0; for (int k = bgs[g]; k < bgs[g + 1]; k++) { btoff[k] = off; off += bl[k]; } }
  const size_t shd = std::max((size_t)tot * 8, (size_t)gw * 8 + (size_t)gw * 2 + 16);
  const double e300 = np1 * 0.5 / K2, e2000 = (np1 + npb) * 0.5 / K2;
  auto cv = [&](double c, double e, int sb) -> unsigned int { if (c <= -1e8) return 0; double v = ceil((c - e) * KC) + sb; return v < 0 ? 0 : (v > 65535 ? 65535 : (unsigned int)v); };
  const unsigned int cut1vd = cv(cut1, errC1 + e300, sumBC1), cut2vd = cv(cut2, errC2 + e2000, sumBC);

  // ---- rows and spans
  struct Span { int row, i0, len; };
  std::vector<RowD> rows; std::vector<Span> spans; std::vector<std::array<long long, 2>> vpairs;
  if (valmode) {
    FILE* fp = fopen(pairsf.c_str(), "r"); if (!fp) { fprintf(stderr, "cannot open %s\n", pairsf.c_str()); return 1; }
    long long pk, qk, M; char ln[256];
    if (g_fam == 7) { while (fgets(ln, sizeof ln, fp)) { long long u, sv, Mv = 1; int nf = sscanf(ln, "%lld %lld %lld", &u, &sv, &Mv); if (nf < 2) continue; if (nf < 3) Mv = 1;
        if (Mv <= 0 || sv % Mv) { fprintf(stderr, "bad z4b pair\n"); return 1; }
        RowD r; r.q = u; r.M = Mv; r.tri = 0; r.pad = 0; rows.push_back(r); spans.push_back({(int)rows.size() - 1, (int)(sv / Mv), 1}); }
      fclose(fp); fp = nullptr; }
    else if (g_fam == 5) { while (fgets(ln, sizeof ln, fp)) { char* e; long long a = strtoll(ln, &e, 10); while (*e == ' ') e++; __int128 b = 0; while (*e >= '0' && *e <= '9') b = b * 10 + (*e++ - '0');
        if (b <= 0 || llabs(a) >= (1LL << 50)) { fprintf(stderr, "bad z3w pair\n"); return 1; }
        RowD r; r.q = (long long)(unsigned long long)b; r.M = (long long)(b >> 64); r.tri = (int)(unsigned int)(unsigned long long)a; r.pad = (int)(a >> 32); rows.push_back(r); spans.push_back({(int)rows.size() - 1, 0, 1}); }
      fclose(fp); fp = nullptr; }
    if (fp)
    while (fgets(ln, sizeof ln, fp)) if (sscanf(ln, "%lld %lld %lld", &pk, &qk, &M) == 3) {
      if (M <= 0 || pk % M || qk <= 0 || qk >= (g_fam == 4 ? (1LL << 62) : (1LL << 40)) || llabs(pk / M) >= (1LL << 31) || (g_conic && M != 1)) { fprintf(stderr, "bad pair %lld %lld %lld\n", pk, qk, M); return 1; }
      RowD r; r.q = qk; r.M = M; r.tri = 0; r.pad = 0; rows.push_back(r); spans.push_back({(int)rows.size() - 1, (int)(pk / M), 1}); }
    if (fp) fclose(fp);
  } else {
    const bool r64 = exists(spansf + ".rows64");
    std::string rb = slurp(spansf + (r64 ? ".rows64" : ".rows")), sb = slurp(spansf + ".spans");
    const size_t rs = r64 ? 24 : 16; size_t nr = rb.size() / rs, ns = sb.size() / 12; rows.resize(nr); spans.resize(ns);
    for (size_t k = 0; k < nr; k++) { memcpy(&rows[k].q, rb.data() + rs * k, 8); rows[k].pad = 0;
      if (r64) { memcpy(&rows[k].M, rb.data() + 24 * k + 8, 8); memcpy(&rows[k].tri, rb.data() + 24 * k + 16, 4); memcpy(&rows[k].pad, rb.data() + 24 * k + 20, 4); }
      else { int m32; memcpy(&m32, rb.data() + 16 * k + 8, 4); rows[k].M = m32; memcpy(&rows[k].tri, rb.data() + 16 * k + 12, 4); } }
    memcpy(spans.data(), sb.data(), 12 * ns);
  }
  double totcand = 0; for (auto& s : spans) totcand += s.len;

  // ---- resume
  std::string ckf = outdir + "/checkpoint";
  size_t ck_sp = 0; long ck_off = 0; double cand_before = 0; unsigned long long surv_before = 0; double wall_before = 0;
  if (resume && exists(ckf)) { FILE* fp = fopen(ckf.c_str(), "r");
    if (fscanf(fp, "%zu %ld %lf %llu %lf", &ck_sp, &ck_off, &cand_before, &surv_before, &wall_before) != 5) { fprintf(stderr, "bad checkpoint\n"); return 1; } fclose(fp);
    fprintf(stderr, "resume at span %zu off %ld\n", ck_sp, ck_off); }

  // ---- watchdog
  signal(SIGTERM, on_sig); signal(SIGINT, on_sig);
  if (!nowd) {
    char pidbuf[32]; snprintf(pidbuf, sizeof pidbuf, "%d", (int)getpid());
    std::string wdpy = std::string(ROOT) + "/search/gpusieve/watchdog.py";
    g_wd = fork();
    if (g_wd == 0) { execlp("python3", "python3", "-W", "ignore", wdpy.c_str(), "--pid", pidbuf, "--tag", runid.c_str(), "--flag", wdflag.c_str(), "--heartbeat", hb.c_str(), (char*)0); _exit(127); }
    atexit(kill_wd);
    double t0 = now(); while (mtime_age(hb) > 5 && now() - t0 < 15) usleep(200000);
    if (mtime_age(hb) > 5) { fprintf(stderr, "watchdog heartbeat missing, abort\n"); return 3; }
  }

  // ---- device
  CK(cudaSetDeviceFlags(cudaDeviceScheduleBlockingSync));
  CK(cudaMemcpyToSymbol(c_l, L1.data(), np1 * 4)); CK(cudaMemcpyToSymbol(c_soff, soff.data(), np1 * 4)); CK(cudaMemcpyToSymbol(c_poff, poff.data(), np1 * 8));
  CK(cudaMemcpyToSymbol(c_np, &np1, 4));
  { int cub = g_fam == 3, nog = g_fam >= 3, swp = g_fam == 4 || g_fam == 5, wid = g_fam == 5; CK(cudaMemcpyToSymbol(c_wide, &wid, 4)); CK(cudaMemcpyToSymbol(c_cube, &cub, 4)); CK(cudaMemcpyToSymbol(c_nogcd, &nog, 4)); CK(cudaMemcpyToSymbol(c_swap, &swp, 4)); }
  CK(cudaMemcpyToSymbol(c2_l, l2.data(), np2 * 4)); CK(cudaMemcpyToSymbol(c2_off, off2.data(), np2 * 4)); CK(cudaMemcpyToSymbol(c2_ioff, ioff2.data(), np2 * 4)); CK(cudaMemcpyToSymbol(c2_m, m2.data(), np2 * 4));
  CK(cudaMemcpyToSymbol(cd_step, dstep.data(), np1 * 4)); CK(cudaMemcpyToSymbol(cd_m, dm.data(), np1 * 4)); CK(cudaMemcpyToSymbol(cd_ioff, dioff.data(), np1 * 4));
  CK(cudaMemcpyToSymbol(cb_l, bl.data(), npb * 4)); CK(cudaMemcpyToSymbol(cb_woff, bwoff.data(), npb * 4)); CK(cudaMemcpyToSymbol(cb_goff, bgoff.data(), npb * 4));
  CK(cudaMemcpyToSymbol(cb_toff, btoff.data(), npb * 4)); CK(cudaMemcpyToSymbol(cb_D, bD.data(), npb * 4)); CK(cudaMemcpyToSymbol(cb_step, bstep.data(), npb * 4));
  CK(cudaMemcpyToSymbol(cb_m, bm.data(), npb * 4)); CK(cudaMemcpyToSymbol(cb_ioff, bioff.data(), npb * 4)); CK(cudaMemcpyToSymbol(cb_gs, bgs.data(), bgs.size() * 4));
  { int ngb = bgs.size() - 1; CK(cudaMemcpyToSymbol(cb_ng, &ngb, 4)); CK(cudaMemcpyToSymbol(cb_np, &npb, 4)); }
  unsigned long long* dPC; CK(cudaMalloc(&dPC, ptot * 8)); CK(cudaMemcpy(dPC, PC.data(), ptot * 8, cudaMemcpyHostToDevice));
  unsigned short* dWB; CK(cudaMalloc(&dWB, WB.size() * 2)); CK(cudaMemcpy(dWB, WB.data(), WB.size() * 2, cudaMemcpyHostToDevice));
  short* dW2; CK(cudaMalloc(&dW2, W2.size() * 2)); CK(cudaMemcpy(dW2, W2.data(), W2.size() * 2, cudaMemcpyHostToDevice));
  unsigned short* dINV; CK(cudaMalloc(&dINV, INV.size() * 2)); CK(cudaMemcpy(dINV, INV.data(), INV.size() * 2, cudaMemcpyHostToDevice));
  if (g_fam == 7) CK(cudaMemcpyToSymbol(g_INVP, &dINV, sizeof(dINV)));
  unsigned int* dSHT; CK(cudaMalloc(&dSHT, SHT.size() * 4)); CK(cudaMemcpy(dSHT, SHT.data(), SHT.size() * 4, cudaMemcpyHostToDevice));
  { int z4f = g_fam == 7; CK(cudaMemcpyToSymbol(c_z4, &z4f, 4)); }
  if (g_fam == 7) {   // z4b tables: kronecker(x, l) and D^-1 mod l at the INV offsets; dense weights of a_l = +1, -1
    std::vector<signed char> QR(INV.size(), 0); std::vector<unsigned short> DV(INV.size(), 0); std::vector<unsigned short> wp(npb), wm(npb);
    for (int k = 0; k < np2; k++) { const long l = l2[k]; const int o = ioff2[k];
      for (long x = 1; x < l; x++) QR[o + x] = -1; for (long x = 1; x < l; x++) QR[o + (x * x) % l] = 1;
      const long dm = ((g_dconst % l) + l) % l; DV[o] = dm ? (unsigned short)modinv(dm, l) : 0; }
    for (int k = 0; k < npb; k++) { wp[k] = (unsigned short)(wc(I2[lv1 + k], 1) + BCb[k]); wm[k] = (unsigned short)(wc(I2[lv1 + k], -1) + BCb[k]); }
    CK(cudaMemcpyToSymbol(cb_wp, wp.data(), 2 * npb)); CK(cudaMemcpyToSymbol(cb_wm, wm.data(), 2 * npb));
    signed char* dQR; unsigned short* dDV; CK(cudaMalloc(&dQR, QR.size())); CK(cudaMemcpy(dQR, QR.data(), QR.size(), cudaMemcpyHostToDevice));
    CK(cudaMalloc(&dDV, 2 * DV.size())); CK(cudaMemcpy(dDV, DV.data(), 2 * DV.size(), cudaMemcpyHostToDevice));
    CK(cudaMemcpyToSymbol(g_QR, &dQR, sizeof(dQR))); CK(cudaMemcpyToSymbol(g_DV, &dDV, sizeof(dDV))); }
  RowD* dRows; CK(cudaMalloc(&dRows, sizeof(RowD) * std::max<size_t>(1, rows.size()))); CK(cudaMemcpy(dRows, rows.data(), sizeof(RowD) * rows.size(), cudaMemcpyHostToDevice));
  const unsigned int CAP1 = (unsigned int)std::min(4.0e8, batch * (cut1 <= -1e8 ? 1.0 : 0.45) + 1e6);
  const unsigned int CAP2 = (unsigned int)std::min(1.0e8, std::max(1e6, CAP1 * (cut2 <= -1e8 ? 1.0 : 0.1)));
  { double need = CAP1 * 24.0 + CAP2 * 48.0;
    FILE* fm = popen("free -b | awk '/Mem:/{print $7}'", "r"); double freeb = 0; if (fm) { if (fscanf(fm, "%lf", &freeb) != 1) freeb = 0; pclose(fm); }
    if (freeb > 0 && freeb - need < 12e9) { fprintf(stderr, "buffers %.1f GB would leave %.1f GB available (< 12 GB, build-lottery): lower --batch\n", need / 1e9, (freeb - need) / 1e9); return 5; } }
  int *dAp, *dAq, *dBp, *dBq, *dBs, *dBa; SieveOut *dC, *dD; unsigned int* dn; Seg6* dseg;
  CK(cudaMalloc(&dAp, 4ull * CAP1)); CK(cudaMalloc(&dAq, 4ull * CAP1)); CK(cudaMalloc(&dBp, 4ull * CAP1)); CK(cudaMalloc(&dBq, 4ull * CAP1)); CK(cudaMalloc(&dBs, 4ull * CAP1)); CK(cudaMalloc(&dBa, 4ull * CAP1));
  CK(cudaMalloc(&dC, sizeof(SieveOut) * (size_t)CAP2)); CK(cudaMalloc(&dD, sizeof(SieveOut) * (size_t)CAP2));
  CK(cudaMalloc(&dn, 16));
  const int MAXSEG = 1 << 21; CK(cudaMalloc(&dseg, sizeof(Seg6) * MAXSEG));
  CK(cudaFuncSetAttribute(dense12, cudaFuncAttributeMaxDynamicSharedMemorySize, (int)shd));
  const size_t sh2 = (size_t)maxg * 2; CK(cudaFuncSetAttribute(stage2g, cudaFuncAttributeMaxDynamicSharedMemorySize, (int)sh2));
  int nsm = 48; CK(cudaDeviceGetAttribute(&nsm, cudaDevAttrMultiProcessorCount, 0));
  int bps = 1; CK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&bps, stage2g, NT2, sh2)); const int grid2 = nsm * std::max(1, bps);
  cudaEvent_t ev[5]; for (auto& e : ev) CK(cudaEventCreate(&e));
  std::vector<SieveOut> hout(CAP2);

  FILE* fs = fopen((outdir + (valmode ? "/valsums.txt" : "/survivors.txt")).c_str(), resume ? "a" : "w");
  FILE* fx = valmode ? nullptr : fopen((outdir + "/dropped.txt").c_str(), resume ? "a" : "w");
  FILE* fpg = fopen((outdir + "/progress.tsv").c_str(), resume ? "a" : "w");
  if (!resume) fprintf(fpg, "utc\tspan_done\tqk_done\tcand_done\tfrac_done\tdense_cand_per_s\ts2_in_per_s\ts3_in_per_s\tcand_per_s_busy\tcand_per_s_wall\tn1\tn300\tn2000\tn_final\tgpu_C\tpower_W\tmax_C\tpaused_s\tduty\twall_s\n");
  { FILE* fp = fopen((outdir + "/params.txt").c_str(), resume ? "a" : "w");
    fprintf(fp, "%s start%s sieve6 tag=%s family=%s table=%s spans=%s pairs=%s rows=%zu spans=%zu total_cand=%.0f lims=%d,%d,%d np1=%d np2=%d groups=%zu K2=%g "
                "KC=%g sumB300=%d sumB2000=%d err300<=%.4f err2000<=%.4f dense groups=%zu lane cuts %u %u cut1=%g cut2=%g cut3=%g hshell=[%g,%g) dedup=%d duty=%g dauto=%g batch=%g code=%s\n",
            utc().c_str(), resume ? "(resume)" : "", tag.c_str(), family.c_str(), table.c_str(), spansf.c_str(), pairsf.c_str(), rows.size(), spans.size(), totcand, lim1, lim2, lim3, np1, np2,
            groups.size(), K2, KC, sumBC1, sumBC, errC1, errC2, bgs.size() - 1, cut1vd, cut2vd, cut1, cut2, cut3, hlo, hhi, (int)dedup, duty, dauto, batch, code.c_str());
    fclose(fp); }

  size_t si = ck_sp; long so = ck_off;
  double cand_done = cand_before, t_start = now() - wall_before, busy = 0, b1 = 0, b2 = 0, b3 = 0, paused = 0, last_prog = 0, cand_run = 0, t_run0 = now();
  unsigned long long nsurv = surv_before, n1tot = 0, n300tot = 0, n2000tot = 0, ndup = 0, nout_shell = 0, nsing = 0;
  long nb = 0;
  std::vector<Seg6> segs; std::vector<int> vi, vr;
  while (si < spans.size() && !g_stop && !exists(stopf)) {
    while ((exists(wdflag) || (!nowd && mtime_age(hb) > 10)) && !g_stop && !exists(stopf)) { usleep(500000); paused += 0.5; }
    if (g_stop || exists(stopf)) break;
    double tb = now(), bc = 0; unsigned int n1;
    size_t sj = si; long sjo = so;
    CK(cudaMemset(dn, 0, 16));
    CK(cudaEventRecord(ev[0]));
    if (valmode) {   // every listed pair goes straight to the exact stages
      vi.clear(); vr.clear();
      while (sj < spans.size() && vi.size() < CAP1) { vi.push_back(spans[sj].i0); vr.push_back(spans[sj].row); sj++; }
      n1 = vi.size(); bc = n1;
      CK(cudaMemcpy(dAp, vi.data(), 4ull * n1, cudaMemcpyHostToDevice)); CK(cudaMemcpy(dAq, vr.data(), 4ull * n1, cudaMemcpyHostToDevice));
      CK(cudaEventRecord(ev[1]));
    } else {
      segs.clear();
      while (sj < spans.size() && bc < batch && segs.size() < (size_t)MAXSEG) {
        const Span& S = spans[sj]; const int take = (int)std::min<long>(SEGD, S.len - sjo);
        Seg6 g; g.row = S.row; g.i0 = S.i0 + (int)sjo; g.len = take; g.pad = 0; segs.push_back(g); bc += take; sjo += take;
        if (sjo >= S.len) { sj++; sjo = 0; }
      }
      CK(cudaMemcpy(dseg, segs.data(), sizeof(Seg6) * segs.size(), cudaMemcpyHostToDevice));
      dense12<<<(unsigned int)segs.size(), NTD, shd>>>(dseg, dRows, dPC, dWB, dINV, cut1vd, cut2vd, gw, dAp, dAq, dn, CAP1);
      CK(cudaGetLastError()); CK(cudaEventRecord(ev[1]));
      CK(cudaMemcpy(&n1, dn, 4, cudaMemcpyDeviceToHost));
      if (n1 > CAP1) { fprintf(stderr, "stage-1 buffer overflow (%u > %u): lower --batch\n", n1, CAP1); return 4; }
    }
    unsigned int hn[4];
    CK(cudaEventRecord(ev[2]));
    stage2g<<<grid2, NT2, sh2>>>(0, dAp, dAq, nullptr, nullptr, n1, groups[0].first, groups[0].second, cut1i, dW2, dINV, dSHT, dRows, dBp, dBq, dBs, nullptr, dn + 1, CAP1);
    CK(cudaGetLastError());
    CK(cudaMemcpy(&hn[1], dn + 1, 4, cudaMemcpyDeviceToHost)); const unsigned int n300 = std::min(hn[1], CAP1);
    for (size_t g = 1; g < groups.size(); g++) {
      const bool last = g + 1 == groups.size();
      stage2g<<<grid2, NT2, sh2>>>(last ? 2 : 1, dBp, dBq, g == 1 ? dBs : dBa, dBa, n300, groups[g].first, groups[g].second, cut2i, dW2, dINV, dSHT, dRows, nullptr, nullptr, dBs, dC, dn + 2, CAP2);
      CK(cudaGetLastError());
    }
    CK(cudaEventRecord(ev[3]));
    CK(cudaMemcpy(&hn[2], dn + 2, 4, cudaMemcpyDeviceToHost));
    if (hn[2] > CAP2) { fprintf(stderr, "stage-2 buffer overflow (%u > %u): raise --cut2 or lower --batch\n", hn[2], CAP2); return 4; }
    const unsigned int n2000 = hn[2];
    stage3<<<std::max(1u, (n2000 + 255) / 256), 256>>>(dC, n2000, lv2, np2, cut3i, dW2, dINV, dSHT, dRows, dD, dn + 3, CAP2);
    CK(cudaGetLastError()); CK(cudaEventRecord(ev[4]));
    CK(cudaMemcpy(hn, dn, 16, cudaMemcpyDeviceToHost));
    unsigned int nfin = std::min(hn[3], CAP2);
    CK(cudaDeviceSynchronize());
    float ms1, ms2, ms3; CK(cudaEventElapsedTime(&ms1, ev[0], ev[1])); CK(cudaEventElapsedTime(&ms2, ev[2], ev[3])); CK(cudaEventElapsedTime(&ms3, ev[3], ev[4]));
    const double dt = now() - tb; busy += dt; b1 += ms1 * 1e-3; b2 += ms2 * 1e-3; b3 += ms3 * 1e-3;
    CK(cudaMemcpy(hout.data(), dD, sizeof(SieveOut) * nfin, cudaMemcpyDeviceToHost));
    unsigned int kept = 0; const double wall = now() - t_start;
    for (unsigned int j = 0; j < nfin; j++) {
      const SieveOut& r = hout[j]; const RowD& R = rows[r.q];
      if (g_fam == 7) {   // z4b: u = row value, s = inner (> 0), q = D s^2; keep gcd(u, D s) = 1; print p = u, q = D s^2
        const long long u = R.q, sv = R.M * (long long)r.p, D = g_dconst; const __int128 q = (__int128)D * sv * sv;
        if (valmode) { fprintf(fs, "%lld %s %d %d %d\n", u, s128(q).c_str(), r.s300, r.s2000, r.s10000); kept++; continue; }
        if (sv <= 0 || std::__gcd((unsigned long long)u, (unsigned long long)llabs(D)) != 1 || std::__gcd((unsigned long long)u, (unsigned long long)sv) != 1 || q + 16 * (__int128)u == 0) { nsing++; continue; }
        const long double a1 = (long double)D * sv, a2 = -(long double)D * u, a3 = -(long double)D * D * sv * u, b2 = a1 * a1 + 4 * a2, b4 = a1 * a3, b6 = a3 * a3;
        const long double c4 = b2 * b2 - 24 * b4, c6 = -b2 * b2 * b2 + 36 * b2 * b4 - 216 * b6, H = std::max(3 * logl(fabsl(c4)), 2 * logl(fabsl(c6))) / 12;
        if (H < hlo || H >= hhi) { nout_shell++; continue; }
        fprintf(fs, "z4 %lld %s %.4f %.4f %.4f H=%.6Lf wall_s=%.1f Ds=%lld,%lld\n", u, s128(q).c_str(), r.s300 / K2, r.s2000 / K2, r.s10000 / K2, H, wall, D, sv);
        kept++; continue;
      }
      if (g_fam == 5) {   // z3w: a = inner, |b| = (M << 64) | q; normalise to a > 0; weighted-primitive; print p = b, q = a^3 (decimal, 128-bit)
        long long a = (((long long)R.pad << 32) | (long long)(unsigned int)R.tri) + r.p; __int128 b = ((__int128)R.M << 64) | (__int128)(unsigned long long)R.q; if (a < 0) { a = -a; b = -b; }
        if (a >= 5000000000000LL) { nsing++; continue; }   // a^3 must stay inside int128 (|a| < 5e12 covers H < 29)
        const __int128 a3 = (__int128)a * a * a;
        if (valmode) { fprintf(fs, "%s %s %d %d %d\n", s128(b).c_str(), s128(a3).c_str(), r.s300, r.s2000, r.s10000); kept++; continue; }
        if (a == 0 || b == 0 || a3 == 27 * b) { nsing++; continue; }
        { long long g = (long long)(b % a); if (g < 0) g = -g; long long x = a; while (g) { long long t = x % g; x = g; g = t; } bool prim = true;   // x = gcd(a, |b|)
          for (long long p = 2; x > 1 && p * p <= x; p++) if (x % p == 0) { while (x % p == 0) x /= p; if (b % ((__int128)p * p * p) == 0) prim = false; }
          if (x > 1 && b % ((__int128)x * x * x) == 0) prim = false;
          if (!prim) { nsing++; continue; } }
        const long double A = a, Bb = (long double)b, c4 = A * (A * A * A - 24 * Bb), c6 = -(A * A * A * A * A * A - 36 * A * A * A * Bb + 216 * Bb * Bb);
        const long double H = std::max(3 * logl(fabsl(c4)), 2 * logl(fabsl(c6))) / 12;
        if (H < hlo || H >= hhi) { nout_shell++; continue; }
        fprintf(fs, "z3 %s %s %.4f %.4f %.4f H=%.6Lf wall_s=%.1f ab=%lld,%s\n", s128(b).c_str(), s128(a3).c_str(), r.s300 / K2, r.s2000 / K2, r.s10000 / K2, H, wall, a, s128(b).c_str());
        kept++; continue;
      }
      const __int128 pk128 = (__int128)R.M * r.p;
      if (pk128 > (__int128)9000000000000000000LL || pk128 < -(__int128)9000000000000000000LL) { nsing++; continue; }
      const long long pk = (long long)pk128, qk = R.q;
      if (g_fam == 3 || g_fam == 4) {   // z3: a = row value, b = M i; z3b: |b| = row value, a = M i (normalised to a > 0); keep weighted-primitive pairs (no p with p | a and p^3 | b); print p = b, q = a^3
      long long a = R.q, b = pk; if (g_fam == 4) { a = pk; b = R.q; if (a < 0) { a = -a; b = -b; } }
      const __int128 a3 = (__int128)a * a * a;
        if (valmode) { fprintf(fs, "%lld %s %d %d %d\n", b, s128(a3).c_str(), r.s300, r.s2000, r.s10000); kept++; continue; }
        if (b == 0 || a3 == (__int128)27 * b) { nsing++; continue; }
        { unsigned long long g = std::__gcd((unsigned long long)llabs(a), (unsigned long long)llabs(b)); bool prim = true;
          for (unsigned long long p = 2; g > 1 && p * p <= g; p++) if (g % p == 0) { while (g % p == 0) g /= p; if (b % (long long)(p * p * p) == 0) prim = false; }
          if (g > 1 && b % (long long)(g * g * g) == 0) prim = false;
          if (!prim) { nsing++; continue; } }
        const long double A = a, Bb = b, c4 = A * (A * A * A - 24 * Bb), c6 = -(A * A * A * A * A * A - 36 * A * A * A * Bb + 216 * Bb * Bb);
        const long double H = std::max(3 * logl(fabsl(c4)), 2 * logl(fabsl(c6))) / 12;
        if (H < hlo || H >= hhi) { nout_shell++; continue; }
        if (dedup) { const int nsub = __builtin_popcount(R.tri); int tri = 0, c = 0; for (int w = 0; w < 12 && c < nsub; w++) if (b % Wp[w] == 0 && a % Wp[w] != 0) { tri |= 1 << w; c++; }
          if (tri != R.tri) { ndup++; continue; } }
        fprintf(fs, "z3 %lld %s %.4f %.4f %.4f H=%.6Lf wall_s=%.1f ab=%lld,%lld\n", b, s128(a3).c_str(), r.s300 / K2, r.s2000 / K2, r.s10000 / K2, H, wall, a, b);
        kept++; continue;
      }
      if (g_conic) {   // (u,v) = (i, row q) -> (P,Q) = (F1, F2); primitive points only, sign Q > 0 (or Q = 0, P > 0)
        const __int128 u = r.p, v = R.q;
        __int128 P = CFi[0] * u * u + CFi[1] * u * v + CFi[2] * v * v, Q = CFi[3] * u * u + CFi[4] * u * v + CFi[5] * v * v;
        if (Q < 0 || (Q == 0 && P < 0)) { P = -P; Q = -Q; }
      if (valmode) { fprintf(fs, "%s %s %d %d %d\n", s128(P).c_str(), s128(Q).c_str(), r.s300, r.s2000, r.s10000); kept++; continue; }
        if (CFi[1] == 0 && CFi[4] == 0 && u < 0) { ndup++; continue; }   // diagonal forms: (u,v) and (-u,v) are the same point
        if (P == 0 || Q == 0 || gcd128(P, Q) != 1) { nsing++; continue; }
        if ((g_fam == 0 && (P == Q || Q == -8 * P)) || (g_fam == 2 && (P == Q || P == -Q))) { nsing++; continue; }
        const long double H = Hfam128(P, Q);
        if (H < hlo || H >= hhi) { nout_shell++; continue; }
        fprintf(fs, "%s %s %s %.4f %.4f %.4f H=%.6Lf wall_s=%.1f uv=%d,%lld\n", family.c_str(), s128(P).c_str(), s128(Q).c_str(), r.s300 / K2, r.s2000 / K2, r.s10000 / K2, H, wall, r.p, R.q);
        kept++; continue;
      }
      if (valmode) { fprintf(fs, "%lld %lld %d %d %d\n", pk, qk, r.s300, r.s2000, r.s10000); kept++; continue; }
      const long double H = Hfam(pk, qk);
      if (singular(pk, qk)) { nsing++; continue; }
      if (H < hlo || H >= hhi) { nout_shell++; fprintf(fx, "%s %lld %lld %.4f %.4f %.4f H=%.6Lf out_of_shell\n", family.c_str(), pk, qk, r.s300 / K2, r.s2000 / K2, r.s10000 / K2, H); continue; }
      if (dedup) {   // canonical n-subset: the n = popcount(row mask) smallest W-primes dividing pk (they do not divide qk)
        const int nsub = __builtin_popcount(R.tri);
        int tri = 0, c = 0; for (int w = 0; w < 12 && c < nsub; w++) if (pk % Wp[w] == 0 && qk % Wp[w] != 0) { tri |= 1 << w; c++; }
        if (tri != R.tri) { ndup++; fprintf(fx, "%s %lld %lld %.4f %.4f %.4f H=%.6Lf duplicate_progression\n", family.c_str(), pk, qk, r.s300 / K2, r.s2000 / K2, r.s10000 / K2, H); continue; }
      }
      fprintf(fs, "%s %lld %lld %.4f %.4f %.4f H=%.6Lf wall_s=%.1f\n", family.c_str(), pk, qk, r.s300 / K2, r.s2000 / K2, r.s10000 / K2, H, wall);
      kept++;
    }
    fflush(fs); if (fx) fflush(fx);
    n1tot += n1; n300tot += hn[1]; n2000tot += hn[2];
    nsurv += kept; cand_done += bc; cand_run += bc; nb++;
    si = sj; so = sjo;
    const long long qdone = rows[spans[si < spans.size() ? si : spans.size() - 1].row].q;
    if (!valmode) { FILE* fp = fopen((ckf + ".tmp").c_str(), "w"); fprintf(fp, "%zu %ld %.0f %llu %.3f\n", si, so, cand_done, nsurv, now() - t_start); fclose(fp); rename((ckf + ".tmp").c_str(), ckf.c_str()); }
    if (now() - last_prog > 10 || si >= spans.size()) {
      double gc, w, mx; last_thermal(gc, w, mx);
      fprintf(fpg, "%s\t%zu\t%lld\t%.0f\t%.6f\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%llu\t%llu\t%llu\t%llu\t%.0f\t%.2f\t%.1f\t%.1f\t%.2f\t%.1f\n", utc().c_str(), si, qdone, cand_done,
              cand_done / (totcand > 0 ? totcand : 1), cand_run / b1, n1tot / std::max(b2, 1e-9), n2000tot / std::max(b3, 1e-9), cand_run / busy, cand_run / (now() - t_run0),
              n1tot, n300tot, n2000tot, nsurv, gc, w, mx, paused, duty, now() - t_start);
      fflush(fpg); last_prog = now();
    }
    if (dauto > 0) {
      static double tsamp = 0, tadj = 0, ssum = 0; static int sn = 0;
      if (now() - tsamp >= 3) { double gc, w, mx; last_thermal(gc, w, mx); if (mx == mx) { ssum += mx; sn++; } tsamp = now(); }
      if (now() - tadj >= 30 && sn > 0) {
        double m = ssum / sn; double old = duty;
        if (m > dauto) duty = std::max(0.15, duty - 0.05); else if (m < dauto - 2) duty = std::min(1.0, duty + 0.05);
        if (duty != old) { FILE* fp = fopen((outdir + "/params.txt").c_str(), "a"); fprintf(fp, "%s duty-auto: mean %.1f C over 30 s -> duty %.2f\n", utc().c_str(), m, duty); fclose(fp); }
        ssum = 0; sn = 0; tadj = now();
      }
    }
    if (duty < 1.0 && duty > 0) { double sl = dt * (1 - duty) / duty; usleep((useconds_t)(sl * 1e6)); }
    if (maxbatches > 0 && nb >= maxbatches) break;
  }
  fclose(fs); if (fx) fclose(fx); fclose(fpg);
  const bool done = si >= spans.size();
  FILE* fp = fopen((outdir + "/params.txt").c_str(), "a");
  fprintf(fp, "%s %s cand_run=%.0f busy_s=%.3f dense_s=%.3f stage2_s=%.3f stage3_s=%.3f wall_run_s=%.2f wall_total_s=%.2f paused_s=%.1f dense_cand_per_s=%.4g s2_in_per_s=%.4g s3_in_per_s=%.4g "
              "cand_per_s_busy=%.4g cand_per_s_wall=%.4g n1=%llu n300=%llu n2000=%llu final=%llu dropped_out_of_shell=%llu dropped_duplicates=%llu singular=%llu\n",
          utc().c_str(), done ? "DONE" : "STOPPED", cand_run, busy, b1, b2, b3, now() - t_run0, now() - t_start, paused, cand_run / b1, n1tot / std::max(b2, 1e-9), n2000tot / std::max(b3, 1e-9),
          cand_run / busy, cand_run / (now() - t_run0), n1tot, n300tot, n2000tot, nsurv, nout_shell, ndup, nsing);
  fclose(fp);
  fprintf(stderr, "%s cand=%.4g dense %.4g cand/s | stage2 %.4g in/s | stage3 %.4g in/s | n1=%llu n300=%llu n2000=%llu final=%llu (out-of-shell %llu, dup %llu) | busy %.4g/s wall %.4g/s paused %.0fs\n",
          done ? "DONE" : "STOPPED", cand_run, cand_run / b1, n1tot / std::max(b2, 1e-9), n2000tot / std::max(b3, 1e-9), n1tot, n300tot, n2000tot, nsurv, nout_shell, ndup,
          cand_run / busy, cand_run / (now() - t_run0), paused);
  kill_wd();
  return done ? 0 : 10;
}
