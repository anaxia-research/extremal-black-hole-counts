"""Explicit constants for Part II."""
from mpmath import mp, mpf, sqrt, exp, log, pi, besseli, zeta, ceil
mp.dps = 40

def A(N):      # main term for kappa(N): (pi/sqrt N) I_1(2 pi sqrt N)
    N = mpf(N); return pi / sqrt(N) * besseli(1, 2 * pi * sqrt(N))
def AJ(N):     # main term for c(N): (2 pi/sqrt N) I_1(4 pi sqrt N)
    N = mpf(N); return 2 * pi / sqrt(N) * besseli(1, 4 * pi * sqrt(N))
DJ1 = mpf('0.0308')       # |c(N)/AJ(N) - 1| <= 0.0308 for all N >= 1 (Rademacher + Weil)
def Aratio(Np, N):        # upper bound A(N')/A(N), N' < N
    Np = mpf(Np); N = mpf(N); return (N / Np) ** (mpf(3) / 4) * exp(-2 * pi * (sqrt(N) - sqrt(Np)))
def Alow(N):              # lower bound for A(N)
    N = mpf(N); return exp(2 * pi * sqrt(N)) / (2 * N ** (mpf(3) / 4)) * (1 - 3 / (8 * pi * sqrt(N)))
def q_bound(K):           # q(K) = sum_{s<=K} d(s) <= sqrt2 exp(pi sqrt(K/2))
    return sqrt(2) * exp(pi * sqrt(mpf(K) / 2))
def sigma_bound(m):
    m = mpf(m); return m * (1 + log(m)) if m > 1 else mpf(1)
def g_plus(Y, kap):       # sum kappa(n) e^{-2 pi n Y}, Y >= 1/2 (tail beyond len(kap) negligible, bounded)
    Y = mpf(Y); s = mpf(0)
    for n in range(1, len(kap)):
        t = kap[n] * exp(-2 * pi * n * Y); s += t
        if n > 50 and t < mpf(10) ** -60: break
    return s
def circle_delta(N, G2):
    """analytic bound for |kappa(N) - A(N)| / A(N) from the circle method with Q = ceil(2 N^{1/4})."""
    N = mpf(N); Q = int(ceil(2 * N ** (mpf(1) / 4)))
    E = exp(2 * pi * N / Q ** 2) * (exp(pi / 2 * (1 + mpf(1) / Q) ** 2) / (pi * N) + exp(pi * Q ** 2 / 8) + G2)
    return E / Alow(N)
