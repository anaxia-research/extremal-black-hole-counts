"""Exact Fourier coefficients c(N) of J(tau) = j(tau) - 744 = q^{-1} + sum_{N>=1} c(N) q^N.

Computed as E4^3 / Delta with FLINT integer power series (exact integers)."""
import flint


def _sigma3(L):
    s = [0] * (L + 1)
    for d in range(1, L + 1):
        d3 = d ** 3
        for m in range(d, L + 1, d):
            s[m] += d3
    return s


def euler_product(L):
    """Coefficients of prod_{n>=1} (1 - q^n) up to q^L (Euler pentagonal theorem)."""
    e = [0] * (L + 1)
    k = 0
    while True:
        g1 = k * (3 * k - 1) // 2
        if g1 > L:
            break
        sg = -1 if k % 2 else 1
        e[g1] = sg
        g2 = k * (3 * k + 1) // 2
        if g2 <= L:
            e[g2] = sg
        k += 1
    return e


def _inverse(f, n):
    g = flint.fmpz_poly([1])
    m = 1
    while m < n:
        m = min(2 * m, n)
        e = (f.truncate(m) * g).truncate(m)
        g = (g * (2 - e)).truncate(m)
    return g


def jcoeffs(L):
    """Return list c with c[N] = coefficient of q^N in J for 1 <= N <= L (c[0] = 0)."""
    n = L + 2
    s3 = _sigma3(n)
    E4 = flint.fmpz_poly([1] + [240 * s3[m] for m in range(1, n)])
    E43 = ((E4 * E4).truncate(n) * E4).truncate(n)
    P = flint.fmpz_poly(euler_product(n))
    P2 = (P * P).truncate(n)
    P4 = (P2 * P2).truncate(n)
    P8 = (P4 * P4).truncate(n)
    P16 = (P8 * P8).truncate(n)
    P24 = (P16 * P8).truncate(n)
    qj = (E43 * _inverse(P24, n)).truncate(n)   # q * j
    c = [int(qj[i + 1]) for i in range(L + 1)]
    c[0] = 0                                      # remove the constant 744
    return c
