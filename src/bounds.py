"""Explicit analytic bounds used in the proofs (high-precision floating point, mpmath)."""
from mpmath import mp, mpf, sqrt, exp, log, pi, zeta, besseli, cosh
mp.dps = 40

S = zeta(mpf(3) / 2) ** 2 - 1            # sum_{c>=2} tau(c) c^{-3/2}
DELTA_C = 2 * S * sqrt(2)                # |delta_N| <= DELTA_C * exp(-2 pi sqrt N)


def delta(N):
    return DELTA_C * exp(-2 * pi * sqrt(mpf(N)))


def A(N):
    """Leading Rademacher term (2 pi / sqrt N) I_1(4 pi sqrt N)."""
    N = mpf(N)
    return 2 * pi / sqrt(N) * besseli(1, 4 * pi * sqrt(N))


def sigma_bound(m):
    return m * (1 + log(m)) if m > 1 else mpf(1)


def p_bound(k):
    """p(k) < exp(pi sqrt(2k/3))."""
    return exp(pi * sqrt(mpf(2) * k / 3))


def eta_div(M):
    """Bound for sum_{a>=2, a^2<=M} (1/a) c(M/a^2) / c(M)."""
    M = mpf(M)
    if M < 4:
        return mpf(0)
    return (sqrt(M) / 2) * (1 + delta(1)) / (1 - delta(M)) * 2 * sqrt(2) * exp(-2 * pi * sqrt(M))


def _half_ratio(M):
    """Upper bound for c(floor(M/2)) / c(M):  c(floor(M/2)) <= (1+delta_1) A(M/2)  (A increasing),
    A(M/2)/A(M) <= 2^{3/4} exp(-4 pi (1 - 1/sqrt 2) sqrt M)."""
    M = mpf(M)
    return (1 + delta(1)) / (1 - delta(M)) * 2 ** (mpf(3) / 4) * exp(-4 * pi * (1 - 1 / sqrt(2)) * sqrt(M))


def small_part(k, t):
    """beta_k-type bound:  W^-_k(t) / (k c(k(t+1)))  <=  p(k) sigma(k/2) c(floor(M/2)) / (k c(M)),  M = k(t+1)."""
    M = mpf(k) * (t + 1)
    return p_bound(k) * sigma_bound(mpf(k) / 2) * _half_ratio(M) / k


def small_part_same(k, t):
    """beta'_k-type bound:  W^-_k(t) / (k c(k t)),  M = k t."""
    M = mpf(k) * t
    return p_bound(k) * sigma_bound(mpf(k) / 2) * _half_ratio(M) / k
