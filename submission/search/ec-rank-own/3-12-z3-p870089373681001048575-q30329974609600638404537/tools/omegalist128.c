/* omegalist128.c (build-lottery, own code): all integers b <= X (X up to ~1e36) with omega(b) >= K and every prime factor <= B,
 * by depth-first search over the primes <= B with pruning (a node needing c more distinct primes is dead if multiplying by the next c
 * primes exceeds X).  Optional cap on prime exponents (EMAX) and on b >= XLO.  Output: decimal, one per line, sorted ascending.
 * usage: omegalist128 B X K OUT [EMAX] [XLO]    (X, XLO decimal)      build: gcc -O2 -o omegalist128 omegalist128.c */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned __int128 u128;
static int P[4000], np, K, EMAX = 99; static u128 X, XLO; static u128* buf; static size_t n, cap;
static u128 parse(const char* s) { u128 v = 0; while (*s >= '0' && *s <= '9') v = v * 10 + (*s++ - '0'); return v; }
static void put(FILE* f, u128 v) { char b[64]; int k = 63; b[k] = 0; if (!v) b[--k] = '0'; while (v) { b[--k] = '0' + (int)(v % 10); v /= 10; } fprintf(f, "%s\n", b + k); }
static void push(u128 v) { if (n == cap) { cap = cap ? 2 * cap : 1 << 20; buf = realloc(buf, cap * sizeof(u128)); if (!buf) { fprintf(stderr, "oom\n"); exit(1); } } buf[n++] = v; }
static int feasible(u128 v, int i, int c) { for (int t = 0; t < c; t++) { if (i + t >= np || v > X / P[i + t]) return 0; v *= P[i + t]; } return 1; }
static void dfs(int k, u128 v, int w, int last, int e) {
  if (w >= K && v >= XLO) push(v);
  for (int i = k; i < np; i++) { if (v > X / P[i]) break;
    const int same = (i == last); if (same && e >= EMAX) continue;
    const int w2 = w + !same, need = K - w2;
    if (need > 0 && !feasible(v * P[i], i + 1, need)) { if (!same) break; else continue; }
    dfs(i, v * P[i], w2, i, same ? e + 1 : 1); } }
static int cmp(const void* a, const void* b) { u128 x = *(const u128*)a, y = *(const u128*)b; return x < y ? -1 : x > y; }
int main(int c, char** v) {
  if (c < 5) { fprintf(stderr, "usage: omegalist128 B X K OUT [EMAX] [XLO]\n"); return 1; }
  int B = atoi(v[1]); X = parse(v[2]); K = atoi(v[3]); if (c > 5) EMAX = atoi(v[5]); XLO = c > 6 ? parse(v[6]) : 0;
  for (int p = 2; p <= B && np < 4000; p++) { int ok = 1; for (int d = 2; d * d <= p; d++) if (p % d == 0) { ok = 0; break; } if (ok) P[np++] = p; }
  for (int i = 0; i < np; i++) { if ((u128)P[i] > X) break; if (K > 1 && !feasible(P[i], i + 1, K - 1)) break; dfs(i, P[i], 1, i, 1); }
  qsort(buf, n, sizeof(u128), cmp);
  FILE* fo = fopen(v[4], "w"); for (size_t i = 0; i < n; i++) put(fo, buf[i]); fclose(fo);
  fprintf(stderr, "omegalist128 B %d K %d EMAX %d: %zu values\n", B, K, EMAX, n); return 0; }
