"""Rigorous upper bounds for |delta_N| (kappa(N) = A(N)(1+delta_N)) and eta_M (divisor terms), all N, M >= 1.

delta_bar(N) >= sup_{N' >= N} |delta_{N'}|,  eta_bar(M) >= sup_{M' >= M} eta_{M'}.
For N <= NX the values are exact (from exact kappa); beyond, analytic majorants that decrease."""
from bounds2 import *
from superext import kappa_list
from jcoeffs import jcoeffs

NX = 3000
kap = kappa_list(NX + 5)
cJ = jcoeffs(NX // 4 + 5)
G2 = exp(4 * pi) + g_plus(2, kap)          # g(2i) = g(i/8)
GPH = g_plus(mpf(1) / 2, kap)              # g^+(i/2)

def delta_plus(N):
    """smooth decreasing majorant of the circle-method bound, N >= 400 (Q = ceil(2 N^{1/4}))."""
    N = mpf(N)
    E = exp(pi * sqrt(N) / 2) * (exp(mpf('1.25') * pi / 2) / (pi * N)
                                 + exp(pi * (2 * N ** (mpf(1) / 4) + 1) ** 2 / 8) + G2)
    return E / Alow(N)

with mp.workdps(160):
    _dex = [None] + [abs(mpf(kap[N]) / A(N) - 1) for N in range(1, NX + 1)]
_dsup = [None] * (NX + 2)
_dsup[NX + 1] = delta_plus(NX + 1)
for N in range(NX, 0, -1):
    _dsup[N] = max(_dex[N], _dsup[N + 1])

def delta_bar(N):
    N = int(N)
    return _dsup[max(N, 1)] if N <= NX else delta_plus(N)

def _eps_abs(a, N):
    return kap[N] if a % 4 else cJ[N]

_eex = [None, mpf(0)]
for M in range(2, NX + 1):
    s = sum(mpf(_eps_abs(a, M // (a * a))) / a for a in range(2, int(M ** 0.5) + 1) if M % (a * a) == 0)
    _eex.append(s / kap[M])

def eta_plus(M):
    """decreasing majorant for M >= NX: sqrt(M)/2 * max((1+d1) A(M/4), (1+DJ1) AJ(M/16)) / ((1-d_{M}) A(M))."""
    M = mpf(M)
    num = max((1 + _dsup[1]) * A(M / 4), (1 + DJ1) * AJ(M / 16))
    return sqrt(M) / 2 * num / ((1 - _dsup[400]) * A(M))

_esup = [None] * (NX + 2)
_esup[NX + 1] = eta_plus(NX + 1)
for M in range(NX, 0, -1):
    _esup[M] = max(_eex[M], _esup[M + 1])

def eta_bar(M):
    M = int(M)
    return _esup[max(M, 1)] if M <= NX else eta_plus(M)

def selfcheck():
    out = []
    # circle-method bound dominates the exact values on [400, NX]
    with mp.workdps(160):
        assert all(_dex[N] <= circle_delta(N, G2) <= delta_plus(N) for N in range(400, NX + 1))
    # decreasing majorants (log-derivative argument in the text; grid confirmation)
    grid = [NX + 1, 4000, 10 ** 4, 10 ** 5, 10 ** 6, 10 ** 8]
    assert all(delta_plus(grid[i + 1]) < delta_plus(grid[i]) for i in range(len(grid) - 1))
    assert all(eta_plus(grid[i + 1]) < eta_plus(grid[i]) for i in range(len(grid) - 1))
    assert all(_eex[M] <= eta_plus(M) for M in range(2000, NX + 1))
    out.append("delta_bar(1)=%s delta_bar(5)=%s delta_bar(400)=%s delta_plus(3001)=%s" %
               tuple(mp.nstr(v, 6) for v in [delta_bar(1), delta_bar(5), delta_bar(400), delta_plus(3001)]))
    out.append("eta_bar(1)=%s eta_bar(5)=%s eta_bar(20)=%s eta_plus(3001)=%s" %
               tuple(mp.nstr(v, 6) for v in [eta_bar(1), eta_bar(5), eta_bar(20), eta_plus(3001)]))
    out.append("g(2i) <= %s, g^+(i/2) = %s" % (mp.nstr(G2, 10), mp.nstr(GPH, 8)))
    return out

if __name__ == "__main__":
    for line in selfcheck():
        print(line)
