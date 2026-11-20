// lotenum (build-lottery): z6enum generalised to the families z6k, z5, leg and to any progression class (h, n):
//   outer qk with exactly h primes of W, inner pk = M i with M = product of an n-subset of the W-primes not dividing qk
//   (canonical subset = the n smallest such primes dividing pk; overlaps evaluated, survivors deduplicated by lotsieve7);
//   leg: pk > 0 only (a, b > 0; both orders a > b and a < b are covered by t = pk/qk in (0, inf)).
// usage: lotenum OUT FAM QMAX HLO HHI h n [--qstride k] [--qmin Q] [--share]
//   --share (z3b only, M = 1): keep a with gcd(|a|, b) > 1. Kept pairs are a subset of the
//   default spans. Refuses to write more than 2e7 kept a's (a full list is not a faster sieve:
//   lotsieve7 stage 1 walks contiguous a).
//   --vmul L (z3b only, build-dside 20 Nov 2026): L = comma list of primes (e.g. 2,3). Row b keeps only a = 0 mod M_b,
//   M_b = product of the l in L with v_l(b) = 1 (the 2-/3-adic class of helpers/round4/build-dside.md). Rows carry M_b and
//   go to OUT.vrows (not OUT.rows): stock lotsieve7 ignores M for z3b, so only lotsieve7v --vmul may read them.
//   OUT.info gets 'vmul L positions_full F positions_kept K tau K/F'.
// Below: the original z6enum header (same algorithm).
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
#include "sharecut.h"
typedef double LD;
static const int W[12] = {7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47};
static int FAM = 0;   // 0 z6k, 1 z5, 2 leg, 3 z3, 4 z3b, 7 z4b
static long double DCONST = 1;   // z4b: D of the region
static bool CONIC = false; static LD CF[6];   // conic mode: (pk,qk) = (F1(u,v), F2(u,v)); outer v, inner u, g(s) = H(F1(s,1),F2(s,1))/2
static LD g0(LD t);
static LD g(LD t) {
  if (!CONIC) return g0(t);
  LD P = CF[0] * t * t + CF[1] * t + CF[2], Q = CF[3] * t * t + CF[4] * t + CF[5];
  // H(P,Q) = log|Q| + g0(P/Q) for Q != 0; use the larger coordinate as scale for stability
  if (fabs(Q) >= fabs(P)) return (log(fabs(Q)) + g0(P / Q)) / 2;
  LD r = Q / P; LD a1, a2, a3, a4 = 0;   // H(P,Q) = log|P| + H(1, r)
  if (FAM == 0) { a1 = r - 2; a2 = -r; a3 = -r * (r - 1); } else if (FAM == 1) { a1 = r - 1; a2 = -r; a3 = -r * r; } else { a1 = 0; a2 = 1 + r * r; a3 = 0; a4 = r * r; }
  LD b2 = a1 * a1 + 4 * a2, b4 = 2 * a4 + a1 * a3, b6 = a3 * a3, c4 = b2 * b2 - 24 * b4, c6 = -b2 * b2 * b2 + 36 * b2 * b4 - 216 * b6;
  return (log(fabs(P)) + std::max(3 * log(fabs(c4)), 2 * log(fabs(c6))) / 12) / 2;
}
static LD g0(LD t) {   // H(t,1): H(pk,qk) = log qk + g0(pk/qk)
  LD a1, a2, a3, a4 = 0;
  if (FAM == 0) { a1 = 1 - 2 * t; a2 = -t; a3 = -t * (1 - t); }
  else if (FAM == 1) { a1 = 1 - t; a2 = -t; a3 = -t; }
  else if (FAM == 3) { a1 = 1; a2 = 0; a3 = t; }
  else if (FAM == 4) { a1 = t; a2 = 0; a3 = 1; }   // z3b: H(a,b) = log|b|/3 + g(a/|b|^(1/3)), g(t) = H(a = t, b = 1)
  else if (FAM == 7) { a1 = DCONST * t; a2 = -DCONST; a3 = -DCONST * DCONST * t; }   // z4b: [D s, -D u, -D^2 s u, 0, 0], H(u,s) = log(u)/2 + g(s/sqrt(u)), g(t) = H(u = 1, s = t)
  else { a1 = 0; a2 = t * t + 1; a3 = 0; a4 = t * t; }
  LD b2 = a1 * a1 + 4 * a2, b4 = 2 * a4 + a1 * a3, b6 = a3 * a3;
  LD c4 = b2 * b2 - 24 * b4, c6 = -b2 * b2 * b2 + 36 * b2 * b4 - 216 * b6;
  return std::max(3 * log(fabs(c4)), 2 * log(fabs(c6))) / 12;
}
static LD gu(LD u) { return g(sinh(u)); }
int main(int argc, char** argv) {
  if (argc < 8) { fprintf(stderr, "usage: lotenum OUT FAM QMAX HLO HHI h n [--qstride k] [--qmin Q] [--share]\n"); return 1; }
  std::string out = argv[1], fam = argv[2]; long long QMAX = atoll(argv[3]); LD HLO = strtod(argv[4], 0), HHI = strtod(argv[5], 0);
  const int HW = atoi(argv[6]), NW = atoi(argv[7]);
  if (fam.compare(0, 6, "conic:") == 0) {   // conic:BASE:a1,b1,c1,a2,b2,c2 ; shell is the true H, halved internally
    CONIC = true; std::string rest = fam.substr(6), b = rest.substr(0, rest.find(':')), cs = rest.substr(rest.find(':') + 1);
    fam = b; for (int k = 0; k < 6; k++) { CF[k] = strtold(cs.c_str(), 0); size_t c = cs.find(','); cs = c == std::string::npos ? "" : cs.substr(c + 1); }
    HLO /= 2; HHI /= 2; }
  if (fam == "z6k") FAM = 0; else if (fam == "z5") FAM = 1; else if (fam == "leg") FAM = 2; else if (fam == "z3") FAM = 3; else if (fam == "z3b") FAM = 4; else if (fam == "z4b") FAM = 7; else { fprintf(stderr, "family?\n"); return 1; }
  if (HW < -2 || HW > 12 || NW < 0 || NW > 12 || (HW >= 0 && HW + NW > 12) || (HW == -1 && NW != 0 && FAM != 7)) { fprintf(stderr, "bad h n (h = -1: every outer value, n = 0; h = -2: every h with n = max(T - h, 0), T = n)\n"); return 1; }
  bool SHARE = false; long qstride = 1; long long qmin = 1; LD TLO = -1e300, THI = 1e300; int HMIN = 0; long long MFIX = 1; bool M64 = false; std::vector<long long> MLIST; std::string QLIST; long long FORCE = 1; std::vector<long long> PL; std::vector<size_t> PCO(9, 0); bool PCLASS = false; std::vector<long long> VL; bool VMUL = false; LD vm_full = 0;   // --plist: sorted signed inner values (pair mode); --pclass: 8 sub-lists by 2|p, 3|p, 5|p
  for (int i = 8; i < argc; i++) { std::string a = argv[i]; if (a == "--qstride") qstride = atol(argv[++i]); else if (a == "--qmin") qmin = atoll(argv[++i]); else if (a == "--hmin") HMIN = atoi(argv[++i]); else if (a == "--mfix") MFIX = atoll(argv[++i]); else if (a == "--m64") M64 = true; else if (a == "--qlist") QLIST = argv[++i]; else if (a == "--force") FORCE = atoll(argv[++i]); else if (a == "--dconst") DCONST = strtold(argv[++i], 0);
    else if (a == "--plist") { FILE* fp = fopen(argv[++i], "rb"); if (!fp) { fprintf(stderr, "cannot open plist\n"); return 1; }
      long long v; while (fread(&v, 8, 1, fp) == 1) PL.push_back(v); fclose(fp);
      for (size_t k = 1; k < PL.size(); k++) if (PL[k] <= PL[k-1]) { fprintf(stderr, "plist not strictly increasing\n"); return 1; }
      if (PL.size() >= (1ULL << 31)) { fprintf(stderr, "plist too long\n"); return 1; } }
    else if (a == "--pclass") PCLASS = true;
    else if (a == "--mlist") { FILE* fm = fopen(argv[++i], "r"); long long v; if (!fm) { fprintf(stderr, "cannot open mlist\n"); return 1; } while (fscanf(fm, "%lld", &v) == 1) MLIST.push_back(v); fclose(fm); M64 = true; } else if (a == "--share") SHARE = true; else if (a == "--vmul") { VMUL = true; std::string v = argv[++i]; size_t c = 0; while (c < v.size()) { VL.push_back(atoll(v.c_str() + c)); size_t d = v.find(',', c); if (d == std::string::npos) break; c = d + 1; } } else if (a == "--tlo") TLO = strtod(argv[++i], 0); else if (a == "--thi") THI = strtod(argv[++i], 0); }
  if (VMUL) { for (long long l : VL) if (l != 2 && l != 3) { fprintf(stderr, "--vmul: primes must be 2 or 3 (table primes start at 5; M_b must stay coprime to them)\n"); return 1; }
    if (FAM != 4 || CONIC || SHARE || MFIX != 1 || FORCE != 1 || !MLIST.empty() || !PL.empty() || HW != -1) { fprintf(stderr, "--vmul needs z3b with h = -1 and no --share/--mfix/--force/--mlist/--plist\n"); return 1; } }
  if (SHARE && FAM != 4) { fprintf(stderr, "--share is only for z3b here; z3w uses lotenumwd --share\n"); return 1; }
  if (SHARE && (MFIX != 1 || FORCE != 1 || !MLIST.empty() || !PL.empty())) { fprintf(stderr, "--share needs M = 1 (no --mfix, --force, --mlist, --plist)\n"); return 1; }
  if (PCLASS && !PL.empty()) {   // reorder: class c = (2|p) + 2 (3|p) + 4 (5|p), each class sorted; write OUT.P for lotpair
    std::vector<long long> R; R.reserve(PL.size());
    for (int c = 0; c < 8; c++) { PCO[c] = R.size(); for (long long v : PL) { long long a = v < 0 ? -v : v; int cc = (a % 2 == 0) | ((a % 3 == 0) << 1) | ((a % 5 == 0) << 2); if (cc == c) R.push_back(v); } }
    PCO[8] = R.size(); PL.swap(R);
    FILE* fp = fopen((out + ".P").c_str(), "wb"); fwrite(PL.data(), 8, PL.size(), fp); fclose(fp);
    for (int c = 0; c < 8; c++) fprintf(stderr, "pclass %d: %zu values\n", c, PCO[c + 1] - PCO[c]); }
  const LD PADH = 1e-4, U = CONIC ? 20 : (FAM == 3 ? 40 : (!PL.empty() ? 32 : 14));   // pair mode: small outer values need |t| up to ~1e13
  // monotone pieces of g(sinh u) on [-U, U]: grid, sign changes of the difference, extrema refined by ternary search
  const int N = CONIC ? (1 << 20) : (1 << 22); std::vector<LD> ext; ext.push_back(-U);
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
  if (HW == -1) {   // every outer value: stop where the shell ends, v < exp(HHI - min g) (min over the extrema and the grid ends)
    LD gmin = 1e30; for (size_t k = 0; k < ext.size(); k++) gmin = std::min(gmin, gu(ext[k]));
    { const int NG = 1 << 16; for (int k = 0; k <= NG; k++) gmin = std::min(gmin, gu(-U + 2 * U * k / NG)); }
    long long vcap = (long long)ceil(exp((FAM == 4 ? 3 : (FAM == 7 ? 2 : 1)) * (HHI - gmin + 1e-3))) + 1;   // z3b: H = log|b|/3 + g; z4b: log(u)/2 + g
    fprintf(stderr, "min g %.6f -> outer values capped at %lld (requested %lld)\n", (double)gmin, vcap, QMAX);
    if (vcap < QMAX) QMAX = vcap; }
  // outer values
  struct Q { long long q; int mask; };
  std::vector<Q> qs;
  if (!QLIST.empty() && QLIST.size() > 4 && QLIST.compare(QLIST.size() - 4, 4, ".bin") == 0) {   // binary int64 qlist
    FILE* fq = fopen(QLIST.c_str(), "rb"); if (!fq) { fprintf(stderr, "cannot open qlist\n"); return 1; } long long v;
    while (fread(&v, 8, 1, fq) == 1) if (v > 0 && v >= qmin && v <= QMAX) { int mk = 0; for (int j = 0; j < 12; j++) if (v % W[j] == 0) mk |= 1 << j; qs.push_back({v, mk}); } fclose(fq); }
  else if (!QLIST.empty()) {   // outer values from a file (first column; '#' lines skipped)
    FILE* fq = fopen(QLIST.c_str(), "r"); if (!fq) { fprintf(stderr, "cannot open qlist\n"); return 1; }
    char ln[4096]; while (fgets(ln, sizeof ln, fq)) { if (ln[0] == '#') continue; long long v = atoll(ln); if (v > 0 && v >= qmin && v <= QMAX) qs.push_back({v, 0}); }
    fclose(fq); }
  else if (HW == -1) { for (long long v = std::max(1LL, qmin); v <= QMAX; v++) qs.push_back({v, 0}); }
  else for (int S = 0; S < 4096; S++) { if (HW >= 0 && __builtin_popcount(S) != HW) continue; if (HW == -2 && __builtin_popcount(S) < HMIN) continue;
    long long PS = 1; for (int j = 0; j < 12; j++) if (S >> j & 1) PS *= W[j];
    if (PS > QMAX) continue;
    for (long long m = 1; m <= QMAX / PS; m++) { bool ok = true; for (int j = 0; j < 12; j++) if (!(S >> j & 1) && m % W[j] == 0) { ok = false; break; }
      if (ok && PS * m >= qmin) qs.push_back({PS * m, S}); } }
  std::sort(qs.begin(), qs.end(), [](const Q& a, const Q& b) { return a.q < b.q; });
  size_t nq_all = qs.size();
  if (qstride > 1) { std::vector<Q> s2; for (size_t k = 0; k < qs.size(); k += qstride) s2.push_back(qs[k]); qs.swap(s2); }
  if (!qs.empty() && qs.back().q >= ((FAM == 4 || FAM == 7) ? (1LL << 62) : (1LL << 40))) { fprintf(stderr, "qk too large\n"); return 1; }
  FILE* fr = fopen((out + (VMUL ? (M64 ? ".vrows64" : ".vrows") : (M64 ? ".rows64" : ".rows"))).c_str(), "wb"); FILE* fs = fopen((out + ".spans").c_str(), "wb");
  long long nrows = 0, nspans = 0; LD npos = 0, nuniq = 0, nclip = 0; long long maxabsi = 0;
  long long share_before = 0, share_kept = 0; bool share_over = false; const long long SHARE_CAP = 20000000;
  for (auto& Qv : qs) {
    const long long qk = Qv.q; const LD lq = FAM == 4 ? log((LD)qk) / 3 : (FAM == 7 ? log((LD)qk) / 2 : log((LD)qk)), A = HLO - lq - PADH, B = HHI - lq + PADH;
    // t-intervals where g in [A, B]
    std::vector<std::pair<LD, LD>> iv;
    for (size_t k = 0; k + 1 < ext.size(); k++) {
      LD ua = ext[k], ub = ext[k + 1], ga = gu(ua), gb = gu(ub); bool inc = gb >= ga;
      LD lo = std::min(ga, gb), hi = std::max(ga, gb);
      if (hi < A || lo > B) continue;
      auto cross = [&](LD level) { LD x = ua, y = ub; for (int it = 0; it < 64; it++) { LD m = (x + y) / 2; if ((gu(m) < level) == inc) x = m; else y = m; } return (x + y) / 2; };
      LD u1 = lo >= A ? (inc ? ua : ub) : cross(A), u2 = hi <= B ? (inc ? ub : ua) : cross(B);
      if (u1 > u2) std::swap(u1, u2);
      if ((k == 0 && u1 <= -U + 1e-9) || (k + 2 == ext.size() && u2 >= U - 1e-9)) { fprintf(stderr, "shell reaches |u| = U at qk=%lld\n", qk); return 1; }
      iv.push_back({sinh(u1), sinh(u2)});
    }
    std::sort(iv.begin(), iv.end());
    std::vector<std::pair<long long, long long>> pv;   // merged pk intervals
    const LD qscale = FAM == 3 ? (LD)qk * (LD)qk * (LD)qk : (FAM == 4 ? cbrtl((LD)qk) : (FAM == 7 ? sqrtl((LD)qk) : (LD)qk));   // z3: t = b / a^3; z3b: t = a / |b|^(1/3); z4b: t = s / sqrt(u)
    for (auto& x : iv) { LD fa = floor(x.first * qscale) - 1, fb = ceil(x.second * qscale) + 1;
      if (fa < -9.0e18 || fb > 9.0e18) { fprintf(stderr, "inner value beyond int64 at qk=%lld\n", qk); return 1; }
      long long a = (long long)fa, b = (long long)fb;
      if ((FAM == 2 || FAM == 7) && !CONIC) { if (b < 1) continue; if (a < 1) a = 1; }   // leg: a > 0; z4b: s > 0
      if (CONIC && CF[1] == 0 && CF[4] == 0) { if (b < 0) continue; if (a < 0) a = 0; }   // diagonal forms: u >= 0 only
      if (TLO > -1e299) { long long c = (long long)ceil(TLO * qscale); if (b < c) continue; if (a < c) a = c; }   // --tlo: t >= TLO
      if (THI < 1e299) { long long c = (long long)floor(THI * qscale); if (a > c) continue; if (b > c) b = c; }   // --thi: t <= THI
      if (!pv.empty() && a <= pv.back().second + 1) pv.back().second = std::max(pv.back().second, b); else pv.push_back({a, b}); }
    long long FX = 1; bool fskip = false;
    if (FORCE > 1) { long long f = FORCE; for (long long p = 2; f > 1; p++) { if (f % p) continue; long long pe = 1; int e = 0; while (f % p == 0) { f /= p; pe *= p; e++; }
        long long q = qk; int v = 0; while (q % p == 0) { q /= p; v++; } if (v == 0) FX *= pe; else if (v < e) { fskip = true; break; } } }
    if (fskip) continue;
    if (SHARE) { int rc = share_restrict(pv, (unsigned __int128)qk, share_before, share_kept, SHARE_CAP);
      if (rc == 1) return 1; if (rc == 2) { share_over = true; break; } }
    int rest[12], nr = 0; for (int j = 0; j < 12; j++) if (!(Qv.mask >> j & 1)) rest[nr++] = j;
    const int NWr = HW == -2 ? std::max(NW - __builtin_popcount(Qv.mask), 0) : NW;   // h = -2: n depends on the outer value's W-count
    const int NT = MLIST.empty() ? (1 << nr) : (int)MLIST.size();
    for (int T = 0; T < NT; T++) { if (MLIST.empty() && __builtin_popcount(T) != NWr) continue;
      long long VM = 1; if (VMUL) for (long long l : VL) { long long q = qk; int v = 0; while (q % l == 0 && v < 2) { q /= l; v++; } if (v == 1) VM *= l; }   // v counts to 2: v == 1 means exactly l || b
      long long Ml = MFIX * FX * VM; int tri = 0, top = -1;
      if (!MLIST.empty()) Ml = MLIST[T];   // --mlist: one row per listed modulus (tri 0, no subset dedup; the ranker drops repeated curves)
      else for (int u = 0; u < nr; u++) if (T >> u & 1) { Ml *= W[rest[u]]; tri |= 1 << rest[u]; top = u; }
      if (!M64 && Ml > 2000000000LL) { fprintf(stderr, "M = %lld exceeds int32 (use --m64)\n", Ml); return 1; }
      if (Ml >= (1LL << 51)) { fprintf(stderr, "M = %lld exceeds 2^51\n", Ml); return 1; }
      const long long M = Ml;
      // primes of rest not in T and smaller than max(T): pk divisible by one of them belongs to an earlier subset
      std::vector<int> ex; for (int u = 0; u < nr; u++) if (!(T >> u & 1) && u < top) ex.push_back(W[rest[u]]);
      int32_t row = (int32_t)nrows; bool any = false;
      for (auto& P : pv) {
        auto cdiv = [](long long a, long long b) { long long q = a / b; if (q * b != a && ((a < 0) == (b < 0))) q++; return q; };
        auto fdiv = [](long long a, long long b) { long long q = a / b; if (q * b != a && ((a < 0) != (b < 0))) q--; return q; };
        long long i0, i1;
        if (PCLASS && !PL.empty()) {   // pair mode with classes: one span per class coprime to qk's {2,3,5} part
          const int sq = (qk % 2 == 0) | ((qk % 3 == 0) << 1) | ((qk % 5 == 0) << 2);
          for (int c = 0; c < 8; c++) { if (c & sq) continue;
            const long long a0 = std::lower_bound(PL.begin() + PCO[c], PL.begin() + PCO[c + 1], P.first) - PL.begin();
            const long long a1 = (long long)(std::upper_bound(PL.begin() + PCO[c], PL.begin() + PCO[c + 1], P.second) - PL.begin()) - 1;
            if (a1 < a0) continue;
            int32_t rec[3] = {row, (int32_t)a0, (int32_t)(a1 - a0 + 1)}; fwrite(rec, 4, 3, fs); nspans++; any = true; npos += a1 - a0 + 1; nuniq += a1 - a0 + 1; }
          continue; }
        if (!PL.empty()) {   // pair mode: indices j with P.first <= PL[j] <= P.second (spans hold index ranges)
          i0 = std::lower_bound(PL.begin(), PL.end(), P.first) - PL.begin(); i1 = (long long)(std::upper_bound(PL.begin(), PL.end(), P.second) - PL.begin()) - 1;
        } else { i0 = cdiv(P.first, M); i1 = fdiv(P.second, M); if (VMUL) vm_full += (LD)(P.second - P.first + 1); }
        if (i1 < i0) continue;
        if (llabs(i0) > 2000000000LL || llabs(i1) > 2000000000LL) {
          if (FAM != 3) { fprintf(stderr, "i out of int32 range at qk=%lld M=%d\n", qk, M); return 1; }
          const long long c0 = std::max(i0, -2000000000LL), c1 = std::min(i1, 2000000000LL); nclip += (LD)(i1 - i0 + 1) - (c1 >= c0 ? (LD)(c1 - c0 + 1) : 0); if (c1 < c0) continue; i0 = c0; i1 = c1; }
        maxabsi = std::max(maxabsi, std::max(llabs(i0), llabs(i1)));
        // spans of at most 2^30 (int32 len)
        for (long long a0 = i0; a0 <= i1; a0 += (1LL << 30)) {
          long long a1 = std::min(i1, a0 + (1LL << 30) - 1), len = a1 - a0 + 1;
          int32_t rec[3] = {row, (int32_t)a0, (int32_t)len}; fwrite(rec, 4, 3, fs); nspans++; any = true; npos += len; }
        for (int s2 = 0; s2 < (1 << ex.size()); s2++) { long long d = 1; int sg = 1; for (size_t e = 0; e < ex.size(); e++) if (s2 >> e & 1) { d *= ex[e]; sg = -sg; }
          nuniq += sg * (LD)(fdiv(i1, d) - cdiv(i0, d) + 1); }
      }
      if (any) { int64_t q64 = qk; if (M64) { int64_t m64 = M; int32_t tp[2] = {tri, 0}; fwrite(&q64, 8, 1, fr); fwrite(&m64, 8, 1, fr); fwrite(tp, 4, 2, fr); }
        else { int32_t mt[2] = {(int32_t)M, tri}; fwrite(&q64, 8, 1, fr); fwrite(mt, 4, 2, fr); } nrows++; }
    }
  }
  fclose(fr); fclose(fs);
  FILE* fi = fopen((out + ".info").c_str(), "w");
  fprintf(fi, "FAM %s h %d n %d hmin %d mfix %lld force %lld tlo %g thi %g QMAX %lld HLO %.4f HHI %.4f qstride %ld qmin %lld\nouter_values_all %zu outer_values_used %zu rows %lld spans %lld\n"
              "positions_evaluated %.6e positions_unique %.6e max_abs_i %lld qk_first %lld qk_last %lld clipped_int32 %.6e\n",
          fam.c_str(), HW, NW, HMIN, MFIX, FORCE, (double)TLO, (double)THI, QMAX, HLO, HHI, qstride, qmin, nq_all, qs.size(), nrows, nspans, npos, nuniq, maxabsi, qs.empty() ? 0 : qs.front().q, qs.empty() ? 0 : qs.back().q, nclip);
  if (VMUL) { fprintf(fi, "vmul"); for (size_t k = 0; k < VL.size(); k++) fprintf(fi, "%s%lld", k ? "," : " ", VL[k]);
    fprintf(fi, " positions_full %.6e positions_kept %.6e tau %.6f\n", (double)vm_full, (double)npos, vm_full > 0 ? (double)(npos / vm_full) : 0.0); }
  if (SHARE) fprintf(fi, "share 1 positions_before %lld positions_kept %lld frac_dropped %.6f cap_hit %d\n",
                     share_before, share_kept, share_before ? 1 - (double)share_kept / (double)share_before : 0, share_over ? 1 : 0);
  fclose(fi);
  fprintf(stderr, "outer %zu (of %zu) rows %lld spans %lld positions %.4e unique %.4e\n", qs.size(), nq_all, nrows, nspans, npos, nuniq);
  if (SHARE) fprintf(stderr, "share before %lld kept %lld frac_dropped %.4f%s\n", share_before, share_kept,
                     share_before ? 1 - (double)share_kept / (double)share_before : 0, share_over ? " CAP" : "");
  if (share_over) { remove((out + (M64 ? ".rows64" : ".rows")).c_str()); remove((out + ".spans").c_str());
    fprintf(stderr, "share: kept positions exceed %lld; spans removed. Stage 1 walks contiguous a, so this flag does not shrink a GPU run.\n", SHARE_CAP); return 3; }
  return 0;
}
