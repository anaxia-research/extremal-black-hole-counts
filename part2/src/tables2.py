"""Tables for Part II: exact N_k(j), beta_k, light-regime ratios, Cardy comparison."""
from superext import *
from mpmath import mp, mpf, sqrt, cosh, cos, log, pi, sinh
mp.dps = 30
kap = kappa_list(6000); c = jcoeffs(1100); d = dlist(4200); tri = triangular(4200)
print("kappa(1..10):", kap[1:11])
print("d(0..15):", d[:16])
print("TABLE exact N_k(j), 1<=k<=8, 1<=j<=5")
for k in range(1, 9):
    N, w = ns_primaries(k, 5, kap, c, d, tri)
    print(k, [N[j] for j in range(1, 6)])
print("TABLE N_k(1), N_k(2), 1<=k<=12 and w_k(1)")
for k in range(1, 13):
    N, w = ns_primaries(k, 2, kap, c, d, tri)
    print(k, N[1], N[2], w[1])
b = beta_list(24)
print("beta_1..24:", b[1:25])
N10, w10 = ns_primaries(10, 90, kap, c, d, tri)
print("N_10(j)/w_10(j), j=1..5:", [mp.nstr(mpf(N10[j]) / w10[j], 10) for j in range(1, 6)], " j=90:", mp.nstr(mpf(N10[90]) / w10[90], 10))
print("max_{2<=t<=90} w10(t-1)/w10(t) =", mp.nstr(max(mpf(w10[t - 1]) / w10[t] for t in range(2, 91)), 8))
def M0(k, j):
    mu = mpf(k) - mpf(1) / 8; mup = mpf(k) - mpf(9) / 8; nu = mpf(j) + mpf(1) / 8
    sec = cosh(2 * pi * sqrt(mup * nu)) if mup > 0 else cos(2 * pi * sqrt(-mup * nu))
    return (cosh(2 * pi * sqrt(mu * nu)) - sec) / sqrt(nu)
print("Cardy comparison N_k(j)/M0 - 1:")
for k in [1, 2, 3, 9]:
    N, w = ns_primaries(k, 60 if k > 1 else 400, kap, c, d, tri)
    js = [1, 2, 5, 10, 20, 40, 60] + ([100, 200, 400] if k == 1 else [])
    print(" k=%d:" % k, [(j, mp.nstr(N[j] / M0(k, j) - 1, 4)) for j in js])
N1, _ = ns_primaries(1, 400, kap, c, d, tri)
print("k=1 log N_1(j) and Cardy exponent at j=100,200,400:", [(j, mp.nstr(log(N1[j]), 10), mp.nstr(2 * pi * sqrt((1 - mpf(1) / 8) * (j + mpf(1) / 8)), 10)) for j in [100, 200, 400]])

# Ramond partition function H_k(x) = (-1)^k Phi_k at the cusp 1/2 = (-1)^k f_k(-24 - 4096 x prod(1+x^n)^24)  [Witten (3.29)]
from flint import fmpz_poly
def fpoly(k):
    xg = fmpz_poly([1, 0] + [kap[n] for n in range(1, 3 * k + 5)])
    pw = [fmpz_poly([1])]
    for r in range(1, k + 1): pw.append(fmpz_poly([int(v) for v in (pw[-1] * xg).coeffs()[:2 * k + 3]]))
    w = ns_w(k, 0, kap, c, d)
    rem = {t: w[t] for t in range(-k, 1)}
    f = [0] * (k + 1)
    for r in range(k, -1, -1):
        cr = [int(v) for v in pw[r].coeffs()] + [0] * (3 * k + 5)
        a = rem[-r]; f[r] = a
        for t in range(-r, 1): rem[t] -= a * cr[t + r]
    return f
LR = 10
PR = fmpz_poly([1])
for n in range(1, LR + 2): PR = PR * fmpz_poly([1] + [0] * (n - 1) + [1]) ** 24
PR = fmpz_poly([int(v) for v in PR.coeffs()[:LR + 1]])
gR = fmpz_poly([-24]) - 4096 * fmpz_poly([0, 1]) * PR
print("Ramond partition functions H_k = sum h_n x^n (h_0 = beta_k):")
badR = []
for k in range(1, 13):
    f = fpoly(k); H = fmpz_poly([0]); pw = fmpz_poly([1])
    for r in range(k + 1):
        H += f[r] * pw; pw = fmpz_poly([int(v) for v in (pw * gR).coeffs()[:LR + 1]])
    h = [(-1) ** k * int(v) for v in H.coeffs()[:LR + 1]] + [0] * (LR + 1)
    assert h[0] == b[k]
    badR += [(k, n) for n in range(1, LR + 1) if h[n] <= 0 or h[n] % 2]
    print(k, h[:5])
print("h_n (1<=n<=10, k<=12) not positive or not even:", badR)

# main term M_k(j) of Theorem B(i) (k = 1: second term omitted)
def Mk(k, j):
    mu = mpf(k) - mpf(1) / 8; mup = mpf(k) - mpf(9) / 8; nu = mpf(j) + mpf(1) / 8
    sec = cosh(2 * pi * sqrt(mup * nu)) if k >= 2 else 0
    return (cosh(2 * pi * sqrt(mu * nu)) - sec) / sqrt(nu)
from mpmath import exp, log10, fabs
print("Cardy table N_k(j)/M_k(j) - 1 (paper's M_k):")
NSk = {}
for k in [1, 2, 3, 5, 9, 10, 11, 12]:
    NSk[k], _ = ns_primaries(k, 400, kap, c, d, tri)
for k in [1, 2, 3, 5, 9, 10]:
    print(" k=%d:" % k, [(j, mp.nstr(NSk[k][j] / Mk(k, j) - 1, 4)) for j in [1, 2, 5, 10, 20, 40, 60]])
mp.dps = 260
worst = max((fabs(NSk[k][j] / Mk(k, j) - 1) * exp(pi * sqrt((k - mpf(1) / 8) * (j + mpf(1) / 8))), k, j)
            for k in NSk for j in range(1, 401))
print("max over k in", sorted(NSk), ", j<=400 of |N/M-1| e^{pi sqrt(mu nu)}:", mp.nstr(worst[0], 5), "at (k,j) =", worst[1:])
mp.dps = 30
print("plot data log10 N_k(j) at j = 1,4,9,...,400:")
for k in [1, 2, 5, 10]:
    print(" k=%d:" % k, " ".join("(%d,%s)" % (s, mp.nstr(log10(NSk[k][s * s]), 6)) for s in range(1, 21)))
bb = beta_list(80)
print("log10|beta_k| k=1..60 (+/-):", " ".join("%d:%s%s" % (k, "+" if bb[k] > 0 else "-", mp.nstr(log10(abs(bb[k])), 5)) for k in range(1, 61)))
small = [k for k in range(1, 81) if ns_primaries(k, 1, kap, c, d, tri)[0][1] <= abs(bb[k])]
print("k <= 80 with N_k(1) <= |beta_k|:", small)
