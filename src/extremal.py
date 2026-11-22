"""Witten's extremal partition functions Z_k and their Virasoro primary counts.

Z_k = sum_{n >= -k} w_k(n) q^n,   eta * Z_k = q^{1/24} (q^{-k} - q^{1-k} + sum_{n>=1} p_k(n) q^n).

Two independent exact computations:
  * hecke_w : w_k(n) = sum_{m=1}^k p2(k-m) sum_{a | (m,n)} (m/a) c(m n / a^2)   (n >= 1)
  * direct_Z: Z_k as the polynomial in J matching the polar part (independent check).
"""
from math import gcd
from jcoeffs import jcoeffs, euler_product


def p2list(K):
    """p2[s] = number of partitions of s with all parts >= 2."""
    p = [1] + [0] * K
    for part in range(2, K + 1):
        for i in range(part, K + 1):
            p[i] += p[i - part]
    return p


def plist(K):
    p = [1] + [0] * K
    for part in range(1, K + 1):
        for i in range(part, K + 1):
            p[i] += p[i - part]
    return p


def _divisors(n):
    out = []
    d = 1
    while d * d <= n:
        if n % d == 0:
            out.append(d)
            if d * d != n:
                out.append(n // d)
        d += 1
    return out


def hecke_w(k, T, c, p2):
    """Dictionary t -> w_k(t) for -k <= t <= T (needs len(c) > k*T)."""
    w = {t: p2[k + t] for t in range(-k, 1)}
    for t in range(1, T + 1):
        tot = 0
        for m in range(1, k + 1):
            g = gcd(m, t)
            sub = 0
            for a in _divisors(g):
                sub += (m // a) * c[m * t // (a * a)]
            tot += p2[k - m] * sub
        w[t] = tot
    return w


def pentagonal(M):
    """List of (j, eps_j) with prod (1-q^n) = sum eps_j q^j, j <= M."""
    e = euler_product(M)
    return [(j, e[j]) for j in range(M + 1) if e[j]]


def primaries(k, T, c, p2, pent):
    """p_k(n) for 1 <= n <= T via p_k(n) = sum_{0<=j<n} eps_j (w_k(n-j) - p2(k+n-j))."""
    w = hecke_w(k, T, c, p2)
    pr = {}
    for n in range(1, T + 1):
        pr[n] = sum(s * (w[n - j] - p2[k + n - j]) for j, s in pent if j < n)
    return pr, w


def direct_Z(k, T, c):
    """Independent computation of w_k(t), -k <= t <= T, by matching the polar part with powers of J."""
    L = T + k + 1
    # J as list indexed by exponent offset 1: J[i] = coeff of q^{i-1}
    Jq = [1, 0] + [c[n] for n in range(1, L)]
    Jq = Jq[:L + 1]

    def mul(a, b):
        out = [0] * (L + 1)
        for i, x in enumerate(a):
            if x:
                for j in range(L + 1 - i):
                    out[i + j] += x * b[j]
        return out
    # powers: q^m J^m, index i -> q^{i-m}
    powers = [[1] + [0] * L]
    for m in range(1, k + 1):
        powers.append(mul(powers[-1], Jq))
    p2 = p2list(k + L)
    target = {t: p2[k + t] for t in range(-k, 1)}
    coef = {}
    resid = dict(target)
    Z = {t: 0 for t in range(-k, T + 1)}
    for m in range(k, 0, -1):
        a = resid[-m]
        coef[m] = a
        for i in range(L + 1):
            t = i - m
            if -k <= t <= T:
                Z[t] += a * powers[m][i]
            if -k <= t <= 0:
                resid[t] -= a * powers[m][i]
    Z[0] += resid[0]
    return Z
