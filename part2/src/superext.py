"""Witten's N=1 extremal NS partition functions and their superprimary counts.

Variable q = e^{2 pi i sigma}, sigma = tau/2, so q = x = e^{pi i tau}.
  g(sigma) = K(2 sigma) = q^{-1} prod_{n>=1} (1+q^{2n-1})^24 - 24 = q^{-1} + sum kappa(N) q^N   (= T_4A)
  Faber polynomials g_m = q^{-m} + O(q); Hecke (replication) formula
      kappa_m(n) = sum_{a | (m,n)} (m/a) eps_a(mn/a^2),
      eps_a = kappa (a odd), (-1)^{N+1} kappa(N) (a = 2 mod 4, T_2B), c(N) (4 | a, J).
  Phi_k = F_{k}(2 sigma) = sum_{m=1}^k d(k-m) g_m + d(k),  d(s) = [q^s] prod (1+q^{2n+1})/(1-q^{2n+2}).
  psi = prod (1-q^{2n})/(1+q^{2n-1}) = sum_{n>=0} (-1)^{n(n+1)/2} q^{n(n+1)/2}.
  Phi_k psi = q^{-k}(1-q) + sum_{j>=1} N_k(j) q^j.
"""
from math import gcd
from flint import fmpz_poly
from jcoeffs import jcoeffs


def _trunc(p, L):
    return fmpz_poly([int(v) for v in p.coeffs()[:L + 1]])


def kappa_list(L):
    """kappa[N] for 0 <= N <= L (kappa[0] = 0), from
    q^{-1} prod (1+q^{2n-1})^24 = q^{-1} [E(q^2)^2 / (E(q) E(q^4))]^24,  E = prod (1-q^n)."""
    from jcoeffs import euler_product, _inverse
    n = L + 3
    e = euler_product(n)
    E1 = fmpz_poly(e)
    E2 = fmpz_poly([e[i // 2] if i % 2 == 0 else 0 for i in range(n)])
    E4 = fmpz_poly([e[i // 4] if i % 4 == 0 else 0 for i in range(n)])
    B = ((E2 * E2).truncate(n) * _inverse((E1 * E4).truncate(n), n)).truncate(n)
    A = B.pow_trunc(24, n)
    kap = [int(A[N + 1]) for N in range(L + 1)]
    kap[0] -= 24
    assert kap[0] == 0
    return kap


def dlist(S):
    """d[s] = [q^s] prod_{n>=1} (1+q^{2n+1})/(1-q^{2n+2}), 0 <= s <= S."""
    d = [0] * (S + 1)
    d[0] = 1
    for n in range(1, S + 1):
        e = 2 * n + 1
        if e <= S:  # multiply by (1+q^e)
            for i in range(S, e - 1, -1):
                d[i] += d[i - e]
        e = 2 * n + 2
        if e <= S:  # divide by (1-q^e)
            for i in range(e, S + 1):
                d[i] += d[i - e]
    return d


def triangular(M):
    """[(T, sign)] with psi = sum sign q^T, T <= M."""
    out, n = [], 0
    while n * (n + 1) // 2 <= M:
        T = n * (n + 1) // 2
        out.append((T, (-1) ** (T % 2)))
        n += 1
    return out


def _divisors(n):
    out, d = [], 1
    while d * d <= n:
        if n % d == 0:
            out.append(d)
            if d * d != n:
                out.append(n // d)
        d += 1
    return out


def eps(a, N, kap, c):
    if a % 2:
        return kap[N]
    if a % 4 == 2:
        return kap[N] if N % 2 else -kap[N]
    return c[N]


def kappa_m(m, n, kap, c):
    tot = 0
    for a in _divisors(gcd(m, n)):
        tot += (m // a) * eps(a, m * n // (a * a), kap, c)
    return tot


def ns_w(k, T, kap, c, d):
    """w[t] = [q^t] Phi_k for -k <= t <= T (needs kap, c up to k*T)."""
    w = {t: d[k + t] for t in range(-k, 1)}
    for t in range(1, T + 1):
        w[t] = sum(d[k - m] * kappa_m(m, t, kap, c) for m in range(1, k + 1))
    return w


def ns_primaries(k, T, kap, c, d, tri):
    """N_k(j), 1 <= j <= T, via N_k(j) = sum_i eps_i (w(j - T_i) - d(k + j - T_i))."""
    w = ns_w(k, T, kap, c, d)
    N = {}
    for j in range(1, T + 1):
        N[j] = sum(s * (w[j - Ti] - d[k + j - Ti]) for Ti, s in tri if Ti < j)
    return N, w


def beta_list(S):
    """Witten's index beta_{k} = (-1)^k F_k(1), 1 <= k <= S, from
    sum_k F_k(1) p^k = D(p) (1 - 24 sum_{n>=1} (2n-1) p^{2n-1}/(1+p^{2n-1}))."""
    d = dlist(S)
    L = [0] * (S + 1)
    L[0] = 1
    e = 1
    while e <= S:
        j = 1
        while j * e <= S:
            L[j * e] -= 24 * e * (-1) ** (j + 1)
            j += 1
        e += 2
    B = [sum(d[i] * L[s - i] for i in range(s + 1)) for s in range(S + 1)]
    return [(-1) ** k * B[k] for k in range(S + 1)]
