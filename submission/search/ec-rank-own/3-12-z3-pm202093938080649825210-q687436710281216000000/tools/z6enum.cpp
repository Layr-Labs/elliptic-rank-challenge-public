// z6enum: enumerate the like-for-like region (FILTER-ACCEPTANCE sec. 6c) as rows and spans for sieve6.
//   Kubert Z/6, (pk,qk) = (p, p+q), integral model [qk-2pk, -pk qk, -pk qk (qk-pk), 0, 0].
//   outer qk in [1, QMAX] with exactly six primes of W = {7,...,47} (any exponents; other primes free);
//   inner pk = M i, M = product of a triple T of the six W-primes not dividing qk (all 20 triples; a pk divisible by
//   >= 4 of them lies in several progressions: the evaluated count includes these overlaps, survivors are deduplicated by
//   assigning pk to the triple of its three smallest such primes);
//   pk of both signs, H(pk,qk) = (1/12) log max(|c4|^3, c6^2) of the integral model in [HLO, HHI) (pk intervals padded
//   by PADH in H and 1 unit at each end; the exact height test runs on the survivors).
// Output: OUT.rows (int64 qk, int32 M, int32 tri = bitmask of T over W), OUT.spans (int32 row, int32 i0, int32 len),
//   OUT.info (counts).  Rows are in ascending qk, then triple order; spans in ascending pk.
// usage: z6enum OUT QMAX HLO HHI [--qstride k] [--qmin Q]
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cmath>
#include <vector>
#include <algorithm>
#include <string>
typedef double LD;
static const int W[12] = {7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47};
static LD g(LD t) {   // H(t,1): H(pk,qk) = log qk + g(pk/qk)
  LD a1 = 1 - 2 * t, a2 = -t, a3 = -t * (1 - t), b2 = a1 * a1 + 4 * a2, b4 = a1 * a3, b6 = a3 * a3;
  LD c4 = b2 * b2 - 24 * b4, c6 = -b2 * b2 * b2 + 36 * b2 * b4 - 216 * b6;
  return std::max(3 * log(fabs(c4)), 2 * log(fabs(c6))) / 12;
}
static LD gu(LD u) { return g(sinh(u)); }
int main(int argc, char** argv) {
  if (argc < 5) { fprintf(stderr, "usage: z6enum OUT QMAX HLO HHI [--qstride k] [--qmin Q]\n"); return 1; }
  std::string out = argv[1]; long long QMAX = atoll(argv[2]); LD HLO = strtod(argv[3], 0), HHI = strtod(argv[4], 0);
  long qstride = 1; long long qmin = 1;
  for (int i = 5; i < argc; i++) { std::string a = argv[i]; if (a == "--qstride") qstride = atol(argv[++i]); else if (a == "--qmin") qmin = atoll(argv[++i]); }
  const LD PADH = 1e-4, U = 14;
  // monotone pieces of g(sinh u) on [-U, U]: grid, sign changes of the difference, extrema refined by ternary search
  const int N = 1 << 22; std::vector<LD> ext; ext.push_back(-U);
  { LD du = 2 * U / N, prev = gu(-U), cur = gu(-U + du); int sprev = cur > prev ? 1 : -1; prev = cur;
    for (int k = 2; k <= N; k++) { LD x = -U + k * du; cur = gu(x); int s = cur > prev ? 1 : (cur < prev ? -1 : sprev);
      if (s != sprev) {   // extremum in [x-2du, x]
        LD lo = x - 2 * du, hi = x; bool mx = sprev > 0;
        for (int it = 0; it < 200; it++) { LD m1 = lo + (hi - lo) / 3, m2 = hi - (hi - lo) / 3; if ((gu(m1) < gu(m2)) == mx) lo = m1; else hi = m2; }
        ext.push_back((lo + hi) / 2); }
      sprev = s; prev = cur; } }
  ext.push_back(U);
  fprintf(stderr, "g: %zu monotone pieces on u in [-%g,%g]; g(+-U) = %.3f %.3f\n", ext.size() - 1, (double)U, (double)U, gu(-U), gu(U));
  for (size_t k = 0; k < ext.size(); k++) fprintf(stderr, "  break u=%.6f t=%.6f g=%.6f\n", ext[k], sinh(ext[k]), gu(ext[k]));
  // outer values
  struct Q { long long q; int mask; };
  std::vector<Q> qs;
  for (int S = 0; S < 4096; S++) { if (__builtin_popcount(S) != 6) continue;
    long long PS = 1; for (int j = 0; j < 12; j++) if (S >> j & 1) PS *= W[j];
    if (PS > QMAX) continue;
    for (long long m = 1; m <= QMAX / PS; m++) { bool ok = true; for (int j = 0; j < 12; j++) if (!(S >> j & 1) && m % W[j] == 0) { ok = false; break; }
      if (ok && PS * m >= qmin) qs.push_back({PS * m, S}); } }
  std::sort(qs.begin(), qs.end(), [](const Q& a, const Q& b) { return a.q < b.q; });
  size_t nq_all = qs.size();
  if (qstride > 1) { std::vector<Q> s2; for (size_t k = 0; k < qs.size(); k += qstride) s2.push_back(qs[k]); qs.swap(s2); }
  if (!qs.empty() && qs.back().q >= (1LL << 40)) { fprintf(stderr, "qk >= 2^40 not supported\n"); return 1; }
  FILE* fr = fopen((out + ".rows").c_str(), "wb"); FILE* fs = fopen((out + ".spans").c_str(), "wb");
  long long nrows = 0, nspans = 0; LD npos = 0, nuniq = 0; long long maxabsi = 0;
  for (auto& Qv : qs) {
    const long long qk = Qv.q; const LD lq = log((LD)qk), A = HLO - lq - PADH, B = HHI - lq + PADH;
    // t-intervals where g in [A, B]
    std::vector<std::pair<LD, LD>> iv;
    for (size_t k = 0; k + 1 < ext.size(); k++) {
      LD ua = ext[k], ub = ext[k + 1], ga = gu(ua), gb = gu(ub); bool inc = gb >= ga;
      LD lo = std::min(ga, gb), hi = std::max(ga, gb);
      if (hi < A || lo > B) continue;
      auto cross = [&](LD level) { LD x = ua, y = ub; for (int it = 0; it < 90; it++) { LD m = (x + y) / 2; if ((gu(m) < level) == inc) x = m; else y = m; } return (x + y) / 2; };
      LD u1 = lo >= A ? (inc ? ua : ub) : cross(A), u2 = hi <= B ? (inc ? ub : ua) : cross(B);
      if (u1 > u2) std::swap(u1, u2);
      if ((k == 0 && u1 <= -U + 1e-9) || (k + 2 == ext.size() && u2 >= U - 1e-9)) { fprintf(stderr, "shell reaches |u| = U at qk=%lld\n", qk); return 1; }
      iv.push_back({sinh(u1), sinh(u2)});
    }
    std::sort(iv.begin(), iv.end());
    std::vector<std::pair<long long, long long>> pv;   // merged pk intervals
    for (auto& x : iv) { long long a = (long long)floor(x.first * qk) - 1, b = (long long)ceil(x.second * qk) + 1;
      if (!pv.empty() && a <= pv.back().second + 1) pv.back().second = std::max(pv.back().second, b); else pv.push_back({a, b}); }
    int rest[6], nr = 0; for (int j = 0; j < 12; j++) if (!(Qv.mask >> j & 1)) rest[nr++] = j;
    for (int x = 0; x < nr; x++) for (int y = x + 1; y < nr; y++) for (int z = y + 1; z < nr; z++) {
      const int M = W[rest[x]] * W[rest[y]] * W[rest[z]], tri = (1 << rest[x]) | (1 << rest[y]) | (1 << rest[z]);
      // primes of rest not in T and smaller than max(T): pk divisible by one of them belongs to an earlier triple
      std::vector<int> ex; for (int u = 0; u < nr; u++) if (u != x && u != y && u != z && rest[u] < rest[z]) ex.push_back(W[rest[u]]);
      int32_t row = (int32_t)nrows; bool any = false;
      for (auto& P : pv) {
        auto cdiv = [](long long a, long long b) { long long q = a / b; if (q * b != a && ((a < 0) == (b < 0))) q++; return q; };
        auto fdiv = [](long long a, long long b) { long long q = a / b; if (q * b != a && ((a < 0) != (b < 0))) q--; return q; };
        long long i0 = cdiv(P.first, M), i1 = fdiv(P.second, M);
        if (i1 < i0) continue;
        if (llabs(i0) > 2000000000LL || llabs(i1) > 2000000000LL) { fprintf(stderr, "i out of int32 range\n"); return 1; }
        maxabsi = std::max(maxabsi, std::max(llabs(i0), llabs(i1)));
        long long len = i1 - i0 + 1; int32_t rec[3] = {row, (int32_t)i0, (int32_t)len}; fwrite(rec, 4, 3, fs); nspans++; any = true; npos += len;
        // unique (canonical) count by inclusion-exclusion over ex: i in [i0,i1] with no e in ex dividing M i (i.e. e | i)
        for (int s = 0; s < (1 << ex.size()); s++) { long long d = 1; int sg = 1; for (size_t e = 0; e < ex.size(); e++) if (s >> e & 1) { d *= ex[e]; sg = -sg; }
          nuniq += sg * (LD)(fdiv(i1, d) - cdiv(i0, d) + 1); }
      }
      if (any) { int64_t q64 = qk; int32_t mt[2] = {M, tri}; fwrite(&q64, 8, 1, fr); fwrite(mt, 4, 2, fr); nrows++; }
    }
  }
  fclose(fr); fclose(fs);
  FILE* fi = fopen((out + ".info").c_str(), "w");
  fprintf(fi, "QMAX %lld HLO %.4f HHI %.4f qstride %ld qmin %lld\nouter_values_all %zu outer_values_used %zu rows %lld spans %lld\n"
              "positions_evaluated %.6e positions_unique %.6e max_abs_i %lld qk_first %lld qk_last %lld\n",
          QMAX, HLO, HHI, qstride, qmin, nq_all, qs.size(), nrows, nspans, npos, nuniq, maxabsi, qs.empty() ? 0 : qs.front().q, qs.empty() ? 0 : qs.back().q);
  fclose(fi);
  fprintf(stderr, "outer %zu (of %zu) rows %lld spans %lld positions %.4e unique %.4e\n", qs.size(), nq_all, nrows, nspans, npos, nuniq);
  return 0;
}
