// gpusieve public types: final survivor record (device and host) and the GPU filter hook.
#pragma once
#include <cstdint>
#define SIEVEOUT_HAS_BOUND 1
// sums are integers at scale K (S = s/K); bound = Selmer bound written by the filter hook (-1 if no hook ran;
// 0..250 bound, 253/254 undecided-and-kept)
struct SieveOut { int32_t p, q; int32_t s300, s2000, s10000; int32_t bound; };
struct SieveFilterCtx { char family[16]; float K; int reserved[14]; };
// Filter hook, called once per batch on the final survivors in DEVICE memory (d_rec[0..n)).
// It may compact d_rec in place and must return the new count (synchronously, i.e. after its
// kernels finished).  stream is reserved (default stream).  Default: weak identity in sieve3.cu;
// link an object that defines a strong `gpusieve_filter` to replace it (see README).
extern "C" unsigned int gpusieve_filter(const SieveFilterCtx* ctx, SieveOut* d_rec, unsigned int n, void* stream);
