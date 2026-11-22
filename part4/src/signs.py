"""Sign pattern (-1)^l c(n,l) > 0 for non-polar coefficients of the nine extremal genera, n <= 6;
and prime factors of the denominators of the basis coordinates."""
from cuspcheck import genus
from exist import solve
def factorint(n):
    f, p = [], 2
    while p*p <= n:
        while n % p == 0: f.append(p); n //= p
        p += 1
    if n > 1: f.append(n)
    return f
for m in [1, 2, 3, 4, 5, 7, 8, 11, 13]:
    phi = genus(m, 12)
    bad = [(k//2, l, c) for k, v in phi.items() if k % 2 == 0 and 2 <= k <= 12 for l, c in v.items()
           if 4*m*(k//2) - l*l > 0 and not ((1 if l % 2 == 0 else -1)*c > 0)]
    (sol, labels, _), rank, J = solve(m)
    primes = sorted(set(p for s in sol if s != 0 for p in factorint(int(s.q))))
    print("m=%d sign failures (n<=6): %s  denominator primes: %s" % (m, bad[:5], primes), flush=True)
