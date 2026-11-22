"""Heavy regime of Part II: circle method for P_k = Phi_k * theta_T (weight 1/2 on Gamma_0(4)+).
Farey order Q = ceil(5r/2), r = (nu/mu)^{1/4}, y = 1/Q^2.  Certifies the constants of Section 'Heavy black holes'."""
from consts2 import *
from superext import dlist
mp.dps = 40
AL = mpf(5) / 2
TAU = mpf('1.42')
d = dlist(400)

def qser(f, x, n=400): return sum(f(s) * x ** s for s in range(n))

# ---------- 1. sanity checks of the transformation laws (numerical, at sample points) ----------
def eta(t):
    qq = exp(2j * pi * t); return exp(2j * pi * t / 24) * mp.nprod(lambda n: 1 - qq ** n, [1, mp.inf])
def thetaT(s): return eta(s) * eta(4 * s) / eta(2 * s)
def gfun(s):
    qq = exp(2j * pi * s); return mp.nprod(lambda n: (1 + qq ** (2 * n - 1)) ** 24, [1, mp.inf]) / qq - 24
def H1(t):
    qq = exp(2j * pi * t); return 24 + 4096 * qq * mp.nprod(lambda n: (1 + qq ** n) ** 24, [1, mp.inf])

def sanity():
    mp.dps = 25
    s = mpf('0.137') + mpf('0.61') * 1j
    a = gfun(-1 / (4 * s)) - gfun(s)                              # W4 invariance of g
    b = thetaT(-1 / (4 * s)) - sqrt(-2j * s) * thetaT(s)           # theta_T(W4 s) = (-2 i s)^{1/2} theta_T(s)
    c = gfun(s / (2 * s + 1)) + H1(2 * s + 1)                      # g(V0 s) = -H_1(2s+1)
    lhs = abs(thetaT(s / (2 * s + 1))); rhs = abs(sqrt(2 * s + 1)) / sqrt(2) * abs(eta(s) * eta(s + mpf(1) / 2) / eta(2 * s))
    mp.dps = 40
    return [abs(a), abs(b), abs(c), abs(lhs - rhs)]

# ---------- 2. constants of the a priori bounds ----------
Theta8 = sum(exp(-pi / 4 * (mpf(n * (n + 1)) / 2 + mpf(1) / 8)) for n in range(80))
D_quarter = qser(lambda s: d[s], exp(-pi / 4))
D_4pi = qser(lambda s: d[s], exp(-4 * pi), 60)
# g_m(2i) = e^{4 pi m} + sum_t kappa_m(t) e^{-4 pi t}, kappa_m(t) <= m (1+eta(1)) (1+delta(1)) A(mt),
# A(N) <= e^{2 pi sqrt N}/2 and 2 pi sqrt(mt) <= pi m/2 + 2 pi t, so the tail is <= c_m e^{4 pi m} with
c_tail = lambda m: m * (1 + eta_bar(1)) * (1 + delta_bar(1)) / 2 * exp(pi * m / 2 - 4 * pi * m) * exp(-2 * pi) / (1 - exp(-2 * pi))
corr = max(c_tail(m) for m in range(1, 60))     # decreasing in m
PHI2i_factor = D_4pi * (1 + corr)          # Phi_k(2i) <= PHI2i_factor * e^{4 pi k}
def Tstar_ratio(k):                         # T*_k / e^{4 pi k}
    return Theta8 * (PHI2i_factor + D_quarter * exp(pi * k / 4 - 4 * pi * k)) + exp(-4 * pi * k)
thetaR = 1 / sqrt(2) * mp.nprod(lambda m: (1 + exp(-2 * pi * m)) / (1 - exp(-2 * pi * m)), [1, mp.inf])
g1p = g_plus(1, kap); g1d = sum(n * kap[n] * exp(-2 * pi * n) for n in range(1, 300))
HRb = 4096 * exp(-2 * pi) * mp.nprod(lambda n: (1 + exp(-2 * pi * n)) ** 24, [1, mp.inf])  # |g(V0 s)+24| <= HRb
PhiR_factor = qser(lambda s: d[s], exp(-2 * pi), 60) * (exp(2 * pi) + g1d) / (exp(2 * pi) - g1p - (24 + HRb))
def TR_ratio(k): return thetaR * PhiR_factor * exp(2 * pi * k - 4 * pi * k)

# ---------- 3. the envelope ----------
def Env(mu, r):
    mu = mpf(mu); r = mpf(r)
    s = mu * r ** 2; nu = mu * r ** 4
    L = (1 - 2 * exp(-pi * r ** 2)) / 2
    B = sqrt(nu) * exp(-2 * pi * s * (1 - 1 / AL ** 2)) / L
    E1 = 2 / (pi * nu) * sqrt((AL * r + 2) / 2) * exp(pi * mu / 2 * (1 + 1 / (AL * r)) ** 2) * B
    E2 = 8 * sqrt(AL * r + 1) * exp(pi * mu * (AL * r + 1) ** 2 / 8) * B
    E3 = 4 * sqrt(AL * r + 1) * TAU * exp(4 * pi * (mu + mpf(1) / 8)) * B
    return E1, E2, E3, E1 + E2 + E3

def dlogs(mu, r):
    """upper bounds for the partial log-derivatives of E1, E2, E3 (algebraic parts bounded above, L-part dropped since it helps)."""
    mu = mpf(mu); r = mpf(r); c = 1 - 1 / AL ** 2
    alg_r = AL / (2 * (AL * r + 1)) + 2 / r
    dr = [-2 / r + AL / (2 * (AL * r + 2)) - pi * mu * (1 + 1 / (AL * r)) / (AL * r ** 2) - 4 * pi * mu * r * c,
          alg_r + pi * mu * AL * (AL * r + 1) / 4 - 4 * pi * mu * r * c,
          alg_r - 4 * pi * mu * r * c]
    dm = [-1 / (2 * mu) + pi / 2 * (1 + 1 / (AL * r)) ** 2 - 2 * pi * r ** 2 * c,
          1 / (2 * mu) + pi * (AL * r + 1) ** 2 / 8 - 2 * pi * r ** 2 * c,
          1 / (2 * mu) + 4 * pi - 2 * pi * r ** 2 * c]
    return dr, dm

def Dgain(mu, r):
    mu = mpf(mu); r = mpf(r); nu = mu * r ** 4
    return mpf('0.999') * pi * sqrt(mu / (nu + 1)) - 1 / (2 * nu)

def jmono(mu0, r0, n=20000, span=20):
    worst = None
    for i in range(n):
        a = r0 + mpf(span) * i / n; b = r0 + mpf(span) * (i + 1) / n
        E = Env(mu0, a)[3]
        m = Dgain(mu0, b) - log((1 + E) / (1 - E))
        worst = m if worst is None or m < worst else worst
    rb = r0 + span
    return worst, Env(mu0, rb)[3]

if __name__ == "__main__":
    print("sanity |g(W4 s)-g(s)|, |theta_T(W4 s)-(-2is)^1/2 theta_T(s)|, |g(V0 s)+H1(2s+1)|, ||theta_T(V0 s)| - ...|:",
          [mp.nstr(v, 3) for v in sanity()])
    print("Theta(1/8) = %s, D(e^-pi/4) = %s, Phi_k(2i)/e^{4 pi k} <= %s" % (mp.nstr(Theta8, 10), mp.nstr(D_quarter, 10), mp.nstr(PHI2i_factor, 10)))
    tr = max(Tstar_ratio(k) for k in range(1, 50))
    print("max_k T*_k e^{-4 pi k} = %s  (<= TAU = %s)" % (mp.nstr(tr, 8), TAU)); assert tr <= TAU
    print("|theta_R| <= %s, |g(V0 s)| <= %s, Phi_R factor %s, max_k T_R e^{-4 pi k} = %s" %
          (mp.nstr(thetaR, 8), mp.nstr(24 + HRb, 6), mp.nstr(PhiR_factor, 6), mp.nstr(max(TR_ratio(k) for k in range(1, 50)), 4)))
    assert max(TR_ratio(k) for k in range(1, 50)) < TAU and exp(2 * pi) - g1p > 24 + HRb
    r8 = mpf(8) ** (mpf(1) / 4); r20 = mpf(20) ** (mpf(1) / 4)
    # monotonicity of the envelope: log-derivatives negative on the whole region
    for (mu, r) in [(mpf(7) / 8, r8), (mpf(7) / 8, r20), (mpf(23) / 8, r8), (mpf(79) / 8, r8)]:
        dr, dm = dlogs(mu, r)
        assert max(dr) < 0 and max(dm) < 0, (mu, r, dr, dm)
    print("envelope log-derivatives negative at the corners (and they decrease further in r and mu)")
    for (lab, mu, r) in [("k>=10, j>=8k", mpf(79) / 8, r8), ("k>=3, j>=8k", mpf(23) / 8, r8), ("k<=2, j>=20k", mpf(7) / 8, r20)]:
        E = Env(mu, r)
        print("  corner %-14s E1=%s E2=%s E3=%s E=%s" % (lab, *[mp.nstr(v, 5) for v in E]))
        assert E[3] < 1
    # monotonicity in j
    for (lab, mu0, r0) in [("k>=10", mpf(79) / 8, r8), ("3<=k<=9", mpf(23) / 8, mpf(3)), ("k<=2", mpf(7) / 8, mpf(4))]:
        w, Eend = jmono(mu0, r0)
        print("  j-monotone %-8s from r=%s (nu >= %s mu): min margin %s; E at r0+20: %s" % (lab, mp.nstr(r0, 6), mp.nstr(r0 ** 4, 6), mp.nstr(w, 5), mp.nstr(Eend, 3)))
        assert w > 0
        rb = r0 + 20
        dr, _ = dlogs(mu0, rb)
        assert max(dr) + 2 / rb < 0                      # E r^2 decreases beyond r0+20
        assert 3 * Eend * rb ** 2 < mpf('0.9') * pi       # 3E < 0.9 pi / r^2 <= Dgain there
    # monotonicity in k (heavy): j >= 8(k+1), resp. j >= 40 for k <= 2
    lowk = (mpf(1) / 2 - exp(-pi * sqrt(mpf(8)))) * exp(pi * sqrt(mpf(8)))
    Emax = Env(mpf(23) / 8, r8)[3]
    print("  k-monotone: N_{k+1}/N_k >= %s / %s" % (mp.nstr(lowk, 6), mp.nstr((1 + Emax) / (1 - Emax), 6)))
    assert lowk > (1 + Emax) / (1 - Emax)
    print("HEAVY REGIME OK")
