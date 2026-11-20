// lotenumw (build-lottery, own code): rows and spans for the Z/3 ladder with 128-bit b (lotsieve7 --family z3w).
//   rows = |b| from a decimal list (ascending), inner a of both signs; y^2 + a x y + b y = x^3; H(a,b) = log|b|/3 + g(a/|b|^(1/3)),
//   g(t) = H of (a = t, b = 1); keep a with H in [HLO, HHI) (padded; the sieve applies the exact test).
// Output: OUT.rows64 (int64 lo, int64 hi, int32 tri = 0, int32 pad = 0), OUT.spans (int32 row, int32 a0, int32 len), OUT.info.
// usage: lotenumw OUT BLIST HLO HHI [--offset] [--share]
//   --share: keep a with gcd(|a|, b) > 1. Kept pairs are a subset of the default spans.
//   Refuses to write more than 2e7 kept a's (lotsieve7 stage 1 walks contiguous a).
//   --vmul L (build-dside 20 Nov 2026): L = comma list from {2, 3}. Row b keeps only a = 0 mod M_b, M_b = product of the l in L
//   with v_l(b) = 1. Spans hold i = a / M_b; rows go to OUT.vrows64 with tri = M_b, pad = 0. Stock lotsieve7 reads tri/pad as
//   the offset A0, so only lotsieve7v --vmul may read these rows. Not with --offset or --share.
// usage: lotenumw OUT BLIST HLO HHI [--offset] [--share] [--vmul 2,3]
#include <cstdio>
#include "sharecut.h"
#include <cstdlib>
#include <cstdint>
#include <cmath>
#include <vector>
#include <algorithm>
#include <string>
typedef double LD; typedef __int128 i128;   // double: interval bounds only (the sieve applies the exact H test; spans are padded)
static LD g(LD t) { LD a1 = t, a3 = 1, b2 = a1 * a1, b4 = a1 * a3, b6 = a3 * a3, c4 = b2 * b2 - 24 * b4, c6 = -b2 * b2 * b2 + 36 * b2 * b4 - 216 * b6;
  return std::max(3 * log(fabs(c4)), 2 * log(fabs(c6))) / 12; }
static LD gu(LD u) { return g(sinh(u)); }
int main(int argc, char** argv) {
  if (argc < 5) { fprintf(stderr, "usage: lotenumw OUT BLIST HLO HHI [--offset] [--share] [--vmul 2,3]\n"); return 1; }
  std::string out = argv[1]; LD HLO = strtod(argv[3], 0), HHI = strtod(argv[4], 0); const LD PADH = 1e-4, U = 40;
  bool OFFSET = false, SHARE = false, VMUL = false; std::vector<long long> VL; LD vm_full = 0;
  for (int i = 5; i < argc; i++) { std::string f = argv[i]; if (f == "--offset") OFFSET = true; else if (f == "--share") SHARE = true;
    else if (f == "--vmul" && i + 1 < argc) { VMUL = true; std::string v = argv[++i]; size_t c = 0; while (c < v.size()) { VL.push_back(atoll(v.c_str() + c)); size_t d = v.find(',', c); if (d == std::string::npos) break; c = d + 1; } } else { fprintf(stderr, "unknown flag %s\n", argv[i]); return 1; } }
  if (VMUL) { if (OFFSET || SHARE) { fprintf(stderr, "--vmul cannot be combined with --offset or --share\n"); return 1; }
    for (long long l : VL) if (l != 2 && l != 3) { fprintf(stderr, "--vmul: primes must be 2 or 3 (table primes start at 5)\n"); return 1; } }
  std::vector<i128> B; { FILE* f = fopen(argv[2], "r"); if (!f) { fprintf(stderr, "cannot open list\n"); return 1; } char ln[256];
    while (fgets(ln, sizeof ln, f)) { if (ln[0] < '0' || ln[0] > '9') continue; i128 v = 0; for (char* p = ln; *p >= '0' && *p <= '9'; p++) v = v * 10 + (*p - '0'); if (v > 0) B.push_back(v); } fclose(f); }
  std::vector<LD> ext; ext.push_back(-U);
  { const int N = 1 << 22; LD du = 2 * U / N, prev = gu(-U), cur = gu(-U + du); int sprev = cur > prev ? 1 : -1; prev = cur;
    for (int k = 2; k <= N; k++) { LD x = -U + k * du; cur = gu(x); int s = cur > prev ? 1 : (cur < prev ? -1 : sprev);
      if (s != sprev) { LD lo = x - 2 * du, hi = x; bool mx = sprev > 0;
        for (int it = 0; it < 200; it++) { LD m1 = lo + (hi - lo) / 3, m2 = hi - (hi - lo) / 3; if ((gu(m1) < gu(m2)) == mx) lo = m1; else hi = m2; }
        ext.push_back((lo + hi) / 2); }
      sprev = s; prev = cur; } }
  ext.push_back(U);
  FILE* fr = fopen((out + (VMUL ? ".vrows64" : ".rows64")).c_str(), "wb"); FILE* fs = fopen((out + ".spans").c_str(), "wb");
  long long nrows = 0, nspans = 0, noff = 0; LD npos = 0, nclip = 0; long long maxa = 0;
  long long share_before = 0, share_kept = 0; bool share_over = false; const long long SHARE_CAP = 20000000;
  for (size_t k = 0; k < B.size(); k++) {
    const i128 b = B[k]; const LD lb = log((LD)b) / 3, A = HLO - lb - PADH, Bh = HHI - lb + PADH, sc = cbrt((LD)b);
    std::vector<std::pair<LD, LD>> iv;
    for (size_t e = 0; e + 1 < ext.size(); e++) {
      LD ua = ext[e], ub = ext[e + 1], ga = gu(ua), gb = gu(ub); bool inc = gb >= ga; LD lo = std::min(ga, gb), hi = std::max(ga, gb);
      if (hi < A || lo > Bh) continue;
      auto cross = [&](LD level) { LD x = ua, y = ub; for (int it = 0; it < 64; it++) { LD m = (x + y) / 2; if ((gu(m) < level) == inc) x = m; else y = m; } return (x + y) / 2; };
      LD u1 = lo >= A ? (inc ? ua : ub) : cross(A), u2 = hi <= Bh ? (inc ? ub : ua) : cross(Bh); if (u1 > u2) std::swap(u1, u2);
      if ((e == 0 && u1 <= -U + 1e-9) || (e + 2 == ext.size() && u2 >= U - 1e-9)) { fprintf(stderr, "shell reaches |u| = U\n"); return 1; }
      iv.push_back({sinh(u1), sinh(u2)}); }
    std::sort(iv.begin(), iv.end());
    std::vector<std::pair<long long, long long>> pv;
    for (auto& x : iv) { LD fa = floor(x.first * sc) - 1, fb = ceil(x.second * sc) + 1;
      if (!OFFSET && (fa < -2.0e9 || fb > 2.0e9)) { const LD c0 = std::max(fa, (LD)-2.0e9), c1 = std::min(fb, (LD)2.0e9); nclip += (fb - fa + 1) - (c1 >= c0 ? c1 - c0 + 1 : 0); if (c1 < c0) continue; fa = c0; fb = c1; }   // inner a is int32: clip and count
      if (fa < -4.9e12 || fb > 4.9e12) { fprintf(stderr, "|a| beyond 4.9e12 at row %zu\n", k); return 1; }
      long long a0 = (long long)fa, a1 = (long long)fb; if (!pv.empty() && a0 <= pv.back().second + 1) pv.back().second = std::max(pv.back().second, a1); else pv.push_back({a0, a1}); }
    if (SHARE) { int rc = share_restrict(pv, (unsigned __int128)b, share_before, share_kept, SHARE_CAP);
      if (rc == 1) return 1; if (rc == 2) { share_over = true; break; } }
    int32_t VM = 1;
    if (VMUL) { for (long long l : VL) { i128 q = b; int v = 0; while (q % l == 0 && v < 2) { q /= l; v++; } if (v == 1) VM *= (int32_t)l; }
      auto cdv = [](long long x, long long d) { long long q = x / d; if (q * d != x && x > 0) q++; return q; };
      auto fdv = [](long long x, long long d) { long long q = x / d; if (q * d != x && x < 0) q--; return q; };
      std::vector<std::pair<long long, long long>> pi;
      for (auto& P : pv) { vm_full += (LD)(P.second - P.first + 1); long long i0 = cdv(P.first, VM), i1 = fdv(P.second, VM); if (i1 >= i0) pi.push_back({i0, i1}); }
      pv.swap(pi); }
    bool any = false; const int64_t blo = (int64_t)(uint64_t)b, bhi = (int64_t)(b >> 64);
    for (auto& P : pv) for (long long s0 = P.first; s0 <= P.second; s0 += (1LL << 30)) { long long s1 = std::min(P.second, s0 + (1LL << 30) - 1);
      maxa = std::max(maxa, std::max(llabs(s0), llabs(s1))); npos += s1 - s0 + 1;
      if (OFFSET && (s0 < -2000000000LL || s1 > 2000000000LL)) {   // own row (b, A0 = s0), inner 0 .. len-1
        if (any) { int32_t tp[2] = {0, 0}; fwrite(&blo, 8, 1, fr); fwrite(&bhi, 8, 1, fr); fwrite(tp, 4, 2, fr); nrows++; any = false; }
        int32_t rec[3] = {(int32_t)nrows, 0, (int32_t)(s1 - s0 + 1)}; fwrite(rec, 4, 3, fs); nspans++;
        int32_t tp[2] = {(int32_t)(uint32_t)(uint64_t)s0, (int32_t)(s0 >> 32)}; fwrite(&blo, 8, 1, fr); fwrite(&bhi, 8, 1, fr); fwrite(tp, 4, 2, fr); nrows++; noff++; continue; }
      int32_t rec[3] = {(int32_t)nrows, (int32_t)s0, (int32_t)(s1 - s0 + 1)}; fwrite(rec, 4, 3, fs); nspans++; any = true; }
    if (any) { int32_t tp[2] = {VMUL ? VM : 0, 0}; fwrite(&blo, 8, 1, fr); fwrite(&bhi, 8, 1, fr); fwrite(tp, 4, 2, fr); nrows++; }
  }
  fclose(fr); fclose(fs);
  FILE* fi = fopen((out + ".info").c_str(), "w");
  fprintf(fi, "FAM z3w list %s HLO %.4f HHI %.4f\nvalues %zu rows %lld spans %lld\npositions_evaluated %.6e max_abs_a %lld clipped_abs_a_over_2e9 %.6e offset_rows %lld\n", argv[2], HLO, HHI, B.size(), nrows, nspans, npos, maxa, nclip, noff);
  if (VMUL) { fprintf(fi, "vmul"); for (size_t k = 0; k < VL.size(); k++) fprintf(fi, "%s%lld", k ? "," : " ", VL[k]);
    fprintf(fi, " positions_full %.6e positions_kept %.6e tau %.6f\n", (double)vm_full, (double)npos, vm_full > 0 ? (double)(npos / vm_full) : 0.0); }
  if (SHARE) fprintf(fi, "share 1 positions_before %lld positions_kept %lld frac_dropped %.6f cap_hit %d\n",
                     share_before, share_kept, share_before ? 1 - (double)share_kept / (double)share_before : 0, share_over ? 1 : 0);
  fclose(fi);
  fprintf(stderr, "values %zu rows %lld spans %lld positions %.4e\n", B.size(), nrows, nspans, npos);
  if (SHARE) fprintf(stderr, "share before %lld kept %lld frac_dropped %.4f%s\n", share_before, share_kept,
                     share_before ? 1 - (double)share_kept / (double)share_before : 0, share_over ? " CAP" : "");
  if (share_over) { remove((out + ".rows64").c_str()); remove((out + ".spans").c_str());
    fprintf(stderr, "share: kept positions exceed %lld; spans removed. Stage 1 walks contiguous a, so this flag does not shrink a GPU run.\n", SHARE_CAP); return 3; }
  return 0; }
