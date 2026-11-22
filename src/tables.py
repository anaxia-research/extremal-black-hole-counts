"""Data for the tables and the figure of the paper."""
import pickle, sys
from mpmath import mp, mpf, sqrt, cosh, log, pi
mp.dps = 60
sys.path.insert(0, 'src')
from jcoeffs import jcoeffs
from extremal import p2list, pentagonal, primaries

L = 6000
c = jcoeffs(L); p2 = p2list(L + 50); pent = pentagonal(L + 50)
PR = {k: primaries(k, L // k, c, p2, pent)[0] for k in [1, 2, 3, 5, 10]}

def M(k, n):
    nu = mpf(n) + mpf(1) / 24; mu = k - mpf(1) / 24; mu2 = k - mpf(25) / 24
    x = 4 * pi * sqrt(mu * nu)
    return sqrt(2 / nu) * (cosh(x) - (cosh(4 * pi * sqrt(mu2 * nu)) if k > 1 else 0))

print("relative error p_k(n)/M_k(n) - 1")
for k in [1, 2, 3, 5, 10]:
    row = []
    for n in [1, 2, 5, 10, 50, 100]:
        row.append(mp.nstr(mpf(PR[k][n]) / M(k, n) - 1, 3))
    print(k, row)
print("figure data: (sqrt(n), log10 p_k(n))")
for k in [1, 2, 5, 10]:
    pts = [(float(sqrt(n)), float(log(PR[k][n], 10))) for n in [1, 2, 4, 9, 16, 25, 36, 49, 64, 81, 100, 144, 196, 256, 324, 400] if n in PR[k]]
    print(k, " ".join("(%.3f,%.3f)" % p for p in pts))
