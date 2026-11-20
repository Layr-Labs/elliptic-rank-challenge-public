// vslice (build-dside): switch a running stock z3w sweep to --vmul spans at its checkpoint.
// usage: vslice VIN OUT STOCK CKPT
//   VIN   prefix of a lotenumwd --vmul enum of the SAME list (VIN.vrows64, VIN.spans)
//   STOCK prefix of the stock enum the sweep runs on (STOCK.rows64, STOCK.spans; offset_rows must be 0)
//   CKPT  the sweep's checkpoint 'si so cand surv wall' (written after a STOP)
// Writes OUT.vrows64 / OUT.spans / OUT.info holding exactly the not-yet-sieved pairs that satisfy the rule:
// rows with b > b_r, and in row r = row of stock span si only a >= a_start = i0(si) + so (spans ascend in a within a row).
// Prints the checkpoint line for the resumed run: '0 0 cand surv wall'.
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>
typedef __int128 i128;
static std::string slurp(const std::string& f) { FILE* fp = fopen(f.c_str(), "rb"); if (!fp) { fprintf(stderr, "cannot open %s\n", f.c_str()); exit(1); }
  std::string s; char buf[1 << 16]; size_t n; while ((n = fread(buf, 1, sizeof buf, fp)) > 0) s.append(buf, n); fclose(fp); return s; }
struct R24 { int64_t lo, hi; int32_t tri, pad; };
static i128 bval(const R24& r) { return ((i128)r.hi << 64) | (i128)(uint64_t)r.lo; }
static std::string s128(i128 v) { if (v == 0) return "0"; bool neg = v < 0; if (neg) v = -v; std::string s; while (v) { s += char('0' + (int)(v % 10)); v /= 10; } if (neg) s += '-'; return std::string(s.rbegin(), s.rend()); }
int main(int argc, char** argv) {
  if (argc != 5) { fprintf(stderr, "usage: vslice VIN OUT STOCK CKPT\n"); return 1; }
  std::string vin = argv[1], out = argv[2], st = argv[3];
  size_t si; long so; double cand, wall; unsigned long long surv;
  { FILE* fp = fopen(argv[4], "r"); if (!fp || fscanf(fp, "%zu %ld %lf %llu %lf", &si, &so, &cand, &surv, &wall) != 5) { fprintf(stderr, "bad checkpoint\n"); return 1; } fclose(fp); }
  std::string srb = slurp(st + ".rows64"), ssb = slurp(st + ".spans");
  const size_t nss = ssb.size() / 12, nsr = srb.size() / 24;
  if (si >= nss) { fprintf(stderr, "checkpoint at end of spans: nothing left\n"); return 2; }
  int32_t sp[3]; memcpy(sp, ssb.data() + 12 * si, 12);
  if ((size_t)sp[0] >= nsr) { fprintf(stderr, "bad stock row\n"); return 1; }
  R24 sr; memcpy(&sr, srb.data() + 24 * (size_t)sp[0], 24);
  if (sr.tri != 0 || sr.pad != 0) { fprintf(stderr, "stock row has an offset A0: not supported\n"); return 1; }
  const i128 b_r = bval(sr); const long long a_start = (long long)sp[1] + so;
  std::string vrb = slurp(vin + ".vrows64"), vsb = slurp(vin + ".spans");
  const size_t nvr = vrb.size() / 24, nvs = vsb.size() / 12;
  std::vector<R24> rows(nvr); memcpy(rows.data(), vrb.data(), 24 * nvr);
  size_t k0 = 0; while (k0 < nvr && bval(rows[k0]) < b_r) k0++;
  for (size_t k = 1; k < nvr; k++) if (bval(rows[k]) <= bval(rows[k - 1])) { fprintf(stderr, "vrows not ascending\n"); return 1; }
  FILE* fr = fopen((out + ".vrows64").c_str(), "wb"); FILE* fs = fopen((out + ".spans").c_str(), "wb");
  long double pos = 0; size_t ns = 0, nr = 0; long long trimmed_first = 0;
  std::vector<int32_t> newid(nvr, -1);
  for (size_t k = k0; k < nvr; k++) { newid[k] = (int32_t)nr++; fwrite(&rows[k], 24, 1, fr); }
  for (size_t j = 0; j < nvs; j++) { int32_t s[3]; memcpy(s, vsb.data() + 12 * j, 12); if ((size_t)s[0] < k0) continue;
    long long i0 = s[1], i1 = (long long)s[1] + s[2] - 1;
    if ((size_t)s[0] == k0 && bval(rows[k0]) == b_r) {   // the checkpoint row: keep a = M i >= a_start
      const long long M = rows[k0].tri; long long c = a_start / M; if (c * M < a_start) c++; if (c * M - M >= a_start) c--;   // c = ceil(a_start / M)
      while (c * M < a_start) c++; while ((c - 1) * M >= a_start) c--;
      if (c > i0) { trimmed_first += (c > i1 ? i1 + 1 : c) - i0; i0 = c; } if (i0 > i1) continue; }
    int32_t o[3] = {newid[s[0]], (int32_t)i0, (int32_t)(i1 - i0 + 1)}; fwrite(o, 4, 3, fs); ns++; pos += i1 - i0 + 1; }
  fclose(fr); fclose(fs);
  FILE* fi = fopen((out + ".info").c_str(), "w");
  fprintf(fi, "FAM z3w vslice of %s at stock checkpoint span %zu off %ld: b_r %s a_start %lld\nrows %zu spans %zu positions_kept %.6Le trimmed_in_first_row %lld\n",
          vin.c_str(), si, so, s128(b_r).c_str(), a_start, nr, ns, pos, trimmed_first);
  fclose(fi);
  fprintf(stderr, "vslice: b_r %s a_start %lld -> rows %zu spans %zu positions %.4Le (first vmul row index %zu)\n", s128(b_r).c_str(), a_start, nr, ns, pos, k0);
  printf("0 0 %.0f %llu %.3f\n", cand, surv, wall);
  return 0; }
