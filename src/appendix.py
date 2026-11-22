"""Data for the appendices: exact primary counts and comparison with the main term."""
import sys
from mpmath import mp, mpf, sqrt, cosh, pi
mp.dps = 60
sys.path.insert(0, 'src')
from jcoeffs import jcoeffs
from extremal import p2list, pentagonal, primaries

L = 400
c = jcoeffs(L); p2 = p2list(L + 50); pent = pentagonal(L + 50)

def M(k, n):
    nu = mpf(n) + mpf(1) / 24; mu = k - mpf(1) / 24; mu2 = k - mpf(25) / 24
    x = 4 * pi * sqrt(mu * nu)
    return sqrt(2 / nu) * (cosh(x) - (cosh(4 * pi * sqrt(mu2 * nu)) if k > 1 else 0))

print("J coefficients c(1..8):", [c[i] for i in range(1, 9)])
pr1, w1 = primaries(1, 12, c, p2, pent)
print("k=1: n, w_1(n), p_1(n), p_1(n)/M_1(n)-1")
for n in range(1, 11):
    print(n, w1[n], pr1[n], mp.nstr(mpf(pr1[n]) / M(1, n) - 1, 4))
print("exact table p_k(n), k=1..6, n=1..5")
for k in range(1, 7):
    pr, w = primaries(k, 6, c, p2, pent)
    print(k, [pr[n] for n in range(1, 6)])
print("w_k(n) and p_k(n)/w_k(n) for k=10, n=1..5")
pr, w = primaries(10, 6, c, p2, pent)
print([mp.nstr(mpf(pr[n]) / w[n], 8) for n in range(1, 6)])
