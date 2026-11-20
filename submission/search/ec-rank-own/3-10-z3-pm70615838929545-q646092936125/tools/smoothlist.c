/* smoothlist.c (build-lottery, own code): sorted B-smooth integers n <= X with omega(n) >= K, optional flag filters
 * (flag bits: 1 = 8 | n, 2 = 9 | n, 4 = 5 | n), written as int64.
 *   smoothlist B X K pos    OUT [FLAGS_EQ]          positive n, flags(n) == FLAGS_EQ if given (-1 = any)
 *   smoothlist B X K signed OUT [FLAGS_SUP LO HI]   +-n with flags(n) containing FLAGS_SUP, LO <= value <= HI
 * build: gcc -O2 -o smoothlist smoothlist.c */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>
static int P[400], np, K; static unsigned long long X; static long long *buf; static size_t n, cap; static int feq = -1, fsup = 0; static long long LO = LLONG_MIN, HI = LLONG_MAX; static int sgn;
static int flags(unsigned long long v) { return (v % 8 == 0) | ((v % 9 == 0) << 1) | ((v % 5 == 0) << 2); }
static void push(long long v) { if (n == cap) { cap = cap ? 2 * cap : 1 << 20; buf = realloc(buf, cap * 8); if (!buf) { fprintf(stderr, "oom\n"); exit(1); } } buf[n++] = v; }
static void take(unsigned long long v) {
  int f = flags(v);
  if (!sgn) { if (feq < 0 || f == feq) push((long long)v); return; }
  if ((f & fsup) != fsup) return;
  if ((long long)v >= LO && (long long)v <= HI) push((long long)v);
  if (-(long long)v >= LO && -(long long)v <= HI) push(-(long long)v); }
// prune: a node needing c more distinct primes from index i on is dead if v * P[i] * ... * P[i+c-1] > X
static int feasible(unsigned long long v, int i, int c) { for (int t = 0; t < c; t++) { if (i + t >= np || v > X / P[i + t]) return 0; v *= P[i + t]; } return 1; }
static void dfs(int k, unsigned long long v, int w, int last) {
  if (w >= K) take(v);
  for (int i = k; i < np; i++) { if (v > X / P[i]) break;
    const int w2 = w + (i != last), need = K - w2;
    if (need > 0 && !feasible(v * P[i], i + 1, need)) { if (i != last) break; else continue; }
    dfs(i, v * P[i], w2, i); } }
static int cmp(const void* a, const void* b) { long long x = *(const long long*)a, y = *(const long long*)b; return x < y ? -1 : x > y; }
int main(int c, char** v) {
  if (c < 6) { fprintf(stderr, "usage: smoothlist B X K pos|signed OUT [FLAGS_EQ | FLAGS_SUP LO HI]\n"); return 1; }
  int B = atoi(v[1]); X = strtoull(v[2], 0, 10); K = atoi(v[3]); sgn = !strcmp(v[4], "signed");
  if (!sgn && c > 6) feq = atoi(v[6]);
  if (sgn && c > 6) fsup = atoi(v[6]); if (sgn && c > 8) { LO = atoll(v[7]); HI = atoll(v[8]); }
  for (int p = 2; p <= B; p++) { int ok = 1; for (int d = 2; d * d <= p; d++) if (p % d == 0) { ok = 0; break; } if (ok) P[np++] = p; }
  if (K <= 0) take(1);
  for (int i = 0; i < np; i++) { if (P[i] > X) break; if (K > 1 && !feasible(P[i], i + 1, K - 1)) break; dfs(i, P[i], 1, i); }
  qsort(buf, n, 8, cmp); size_t m = 0; for (size_t i = 0; i < n; i++) if (!m || buf[i] != buf[m - 1]) buf[m++] = buf[i];
  FILE* fo = fopen(v[5], "wb"); fwrite(buf, 8, m, fo); fclose(fo);
  fprintf(stderr, "smoothlist B %d X %llu K %d %s flags %d range [%lld,%lld]: %zu values\n", B, X, K, sgn ? "signed" : "pos", sgn ? fsup : feq, LO, HI, m); return 0; }
