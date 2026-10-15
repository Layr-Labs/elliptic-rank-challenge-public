"""The one prime sieve (fingerprint primes below 1000, minimality primes below 10^6)."""

from __future__ import annotations


def primes_upto(n: int) -> list:
    """Every prime p <= n, ascending."""
    s = bytearray([1]) * (n + 1)
    s[0] = s[1] = 0
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]
