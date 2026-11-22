"""Heavy regime (n >= 8k): circle-method error envelope, positivity, monotonicity in n and k."""
from mpmath import mp, mpf, sqrt, exp, log, pi, cosh, ceil
mp.dps = 40

R0 = mpf(8) ** (mpf(1) / 4)


def Tstar(mu):
    k = mu + mpf(1) / 24
    return mpf('1.003') * exp(4 * pi * k) + mpf('1.002') * exp(pi * k) + 1


def E_env(mu, r):
    """Envelope E1+E2+E3 of the relative error, valid for every integer N in [r, r+1)."""
    c = mpf('0.4999') * sqrt(2)
    E1 = 2 * sqrt(r + 2) / (mpf('0.4999') * pi * sqrt(2) * sqrt(mu) * r ** 2) * exp(-2 * pi * mu * (r ** 2 - (1 + 1 / r) ** 2))
    E2 = 8 * sqrt(r + 1) * sqrt(mu) * r ** 2 / c * exp(-(pi * mu / 2) * (3 * r + 1) * (r - 1))
    E3 = 4 * sqrt(r + 1) * sqrt(mu) * r ** 2 / c * Tstar(mu) * exp(-2 * pi * mu * r ** 2)
    return E1 + E2 + E3, (E1, E2, E3)


mu0 = mpf(23) / 24
E0, parts = E_env(mu0, R0)
print("corner mu=23/24, r=8^(1/4): E1,E2,E3 =", [mp.nstr(p, 6) for p in parts], " total", mp.nstr(E0, 6))
assert E0 < 1

# monotonicity of the envelope in r and mu (sampled check of the analytic claim)
prev = None
for i in range(0, 4001):
    r = R0 + mpf(i) / 200
    e = E_env(mu0, r)[0]
    assert prev is None or e < prev
    prev = e
for j in range(1, 400):
    assert E_env(mu0 + mpf(j) / 40, R0)[0] < E_env(mu0 + mpf(j - 1) / 40, R0)[0]
print("envelope decreasing in r (r0..r0+20) and in mu (23/24..~11): ok")

# Theorem A heavy: p >= (1 - E) M > 0.

# Monotonicity in n: need  D(mu, r) > log((1+E(r))/(1-E(r)))  with
#   D = 1.998 pi sqrt(mu/(nu+1)) - 1/(2 nu),  nu = mu r^4   (lower bound for log M(nu+1) - log M(nu))
def D(mu, r):
    nu = mu * r ** 4
    return mpf('1.998') * pi * sqrt(mu / (nu + 1)) - 1 / (2 * nu)

worst = None
steps = 20000
for i in range(steps):
    r0 = R0 + mpf(i) * 20 / steps
    r1 = R0 + mpf(i + 1) * 20 / steps
    e = E_env(mu0, r0)[0]
    g = D(mu0, r1) - log((1 + e) / (1 - e))
    if worst is None or g < worst[0]:
        worst = (g, r0)
print("n-monotone margin min over r in [r0, r0+20] at mu=23/24:", mp.nstr(worst[0], 6), "at r =", mp.nstr(worst[1], 6))
assert worst[0] > 0
# r >= r0 + 20: E <= E_env(mu0, r0+20) * (tiny), D >= 1.99 pi / (r^2 sqrt(1 + 1/nu)) - ...; ratio check
rr = R0 + 20
print("at r = r0+20: E_env =", mp.nstr(E_env(mu0, rr)[0], 6), " D =", mp.nstr(D(mu0, rr), 6))
# beyond: E_env * r^2 decreasing (checked) and D * r^2 increasing to 1.998 pi / sqrt(mu)... :
for r in [rr, 2 * rr, 10 * rr, 100 * rr]:
    print("   r =", mp.nstr(r, 5), " E*r^2 =", mp.nstr(E_env(mu0, r)[0] * r ** 2, 5), " D*r^2 =", mp.nstr(D(mu0, r) * r ** 2, 5))
# D is increasing in mu and E decreasing in mu, so mu = 23/24 is the worst case.
for mu in [mpf(47) / 24, mpf(71) / 24, mpf(10)]:
    print("   mu =", mp.nstr(mu, 4), " corner margin:", mp.nstr(D(mu, R0) - log((1 + E_env(mu, R0)[0]) / (1 - E_env(mu, R0)[0])), 6))

# Monotonicity in k (both in heavy regime, n >= 8(k+1)): M_{k+1}/M_k >= 0.4999 exp(2 pi sqrt(nu/(mu+1))) >= 0.4999 exp(2 pi sqrt 8)
print("k-monotone: 0.4999 exp(2 pi sqrt(8)) =", mp.nstr(mpf('0.4999') * exp(2 * pi * sqrt(8)), 6), "> (1+E)/(1-E) <=", mp.nstr((1 + E0) / (1 - E0), 6))

# Theorem B(i), second bound: for n >= 64k (r >= 8^{1/2}),  E_env <= exp(-sqrt(mu nu)) = exp(-mu r^2).
r64 = sqrt(mpf(8))
val = E_env(mu0, r64)[0] * exp(mu0 * r64 ** 2)
print("E_env * exp(mu r^2) at mu=23/24, r=sqrt(8):", mp.nstr(val, 6), "(must be < 1; decreasing in r and mu)")
assert val < 1
prev = None
for i in range(2001):
    r = r64 + mpf(i) / 100
    v = E_env(mu0, r)[0] * exp(mu0 * r ** 2)
    assert prev is None or v < prev
    prev = v
print("   decreasing on r in [sqrt8, sqrt8+20]: ok")

# a priori constant: Z_k(2i) <= 1.003 e^{4 pi k}
from mpmath import mpf as _m
def p2(K):
    p = [1] + [0] * K
    for part in range(2, K + 1):
        for i in range(part, K + 1):
            p[i] += p[i - part]
    return p
P2 = p2(400)
def sigma(m):
    return sum(d for d in range(1, m + 1) if m % d == 0)
worst = 0
for k in range(1, 301):
    Z = P2[k] + sum(P2[k - m] * (exp(4 * pi * m) + sigma(m) * (exp(pi * m) + 449)) for m in range(1, k + 1))
    worst = max(worst, Z / exp(4 * pi * k))
print("max_{k<=300} Z_k(2i)-bound / e^{4 pi k} =", mp.nstr(worst, 8), "(<= 1.003; for k>300 the extra terms are < k^3 e^{pi sqrt(2k/3) - 3 pi k})")
assert worst <= mpf('1.003')
