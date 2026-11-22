"""Witten's Ramond index beta_k for N=1 extremal theories: sign pattern for every k and growth.
B(p) = sum_k F_k(1) p^k = D(p) Lambda(p) = (1-p) L(sigma) q^{1/8} / theta_T(sigma) at p = q,
with L = -(1/2 pi i) d/dsigma log(g+24) in M_2(Gamma_0(4)), f = L/theta_T of weight 3/2 and
f(W4 s) = -(-2 i s)^{3/2} f(s).  bhat_k = sum_{i<=k} B_i are the coefficients of f; beta_k = (-1)^k B_k."""
from bounds2 import *
from mpmath import zeta
from superext import beta_list
from mpmath import sinh
mp.dps = 40
ZETA32 = zeta(mpf(3) / 2)

def S1(Y):   # sum_m sigma(m) e^{-2 pi m Y} = sum_n n x^n/(1-x^n), x = e^{-2 pi Y}
    x = exp(-2 * pi * Y); return sum(n * x ** n / (1 - x ** n) for n in range(1, 400))
def Qabs(Y):  # 1/psi = prod (1+x^{2n-1})/(1-x^{2n}) at x = e^{-2 pi Y}
    x = exp(-2 * pi * Y); return mp.nprod(lambda n: (1 + x ** (2 * n - 1)) / (1 - x ** (2 * n)), [1, mp.inf])
# remainder at the cusp infinity, Im u >= 1/8:  |f(u) - e(-u/8)| <= e^{pi Y/4} [ (1+24 S1(Y)) Q(Y) - 1 ],  Y = 1/8
Y8 = mpf(1) / 8
T_inf = exp(pi * Y8 / 4) * ((1 + 24 * S1(Y8)) * Qabs(Y8) - 1)
# cusp 1/2, Im s >= 1/2:  |L_R| <= 2 + 48 sum n x^n/(1-x^n) (x = e^{-4 pi Y}),  |1/theta_R| <= sqrt2 prod (1+e^{-2 pi m})/(1-e^{-2 pi m})
xr = exp(-2 * pi)
T_R = (2 + 48 * sum(n * xr ** n / (1 - xr ** n) for n in range(1, 200))) * sqrt(2) * mp.nprod(lambda m: (1 + exp(-2 * pi * m)) / (1 - exp(-2 * pi * m)), [1, mp.inf])
T = max(T_inf, T_R)

def Ebeta(nu):
    """relative error of bhat_k = -2 sqrt2 sinh(pi sqrt(nu/2)) (1 + theta E), nu = k - 1/8, Q = ceil((128 nu)^{1/4})."""
    nu = mpf(nu); Q = int(ceil((128 * nu) ** (mpf(1) / 4))); ey = exp(2 * pi * nu / Q ** 2)
    R1 = (Q + 1) ** (mpf(3) / 2) * 2 ** (-mpf(3) / 2) * exp(pi / 16 * (1 + mpf(1) / Q) ** 2) * ey / (pi * nu)
    R2 = 2 * Q ** 2 * (ZETA32 - 1) * exp(pi * Q ** 2 / 64) * ey
    R3 = 2 * Q ** 2 * ZETA32 * T * ey
    main_low = sqrt(2) * exp(pi * sqrt(nu / 2)) * (1 - exp(-2 * pi * sqrt(nu / 2)))
    return (R1 + R2 + R3) / main_low

def Eplus(nu):
    """smooth decreasing majorant of Ebeta (Q <= (128 nu)^{1/4} + 1, Q >= (128 nu)^{1/4})."""
    nu = mpf(nu); q0 = (128 * nu) ** (mpf(1) / 4); Q1 = q0 + 1
    ey = exp(2 * pi * nu / q0 ** 2)
    R1 = (Q1 + 1) ** (mpf(3) / 2) * 2 ** (-mpf(3) / 2) * exp(pi / 16 * (1 + 1 / q0) ** 2) * ey / (pi * nu)
    R2 = 2 * Q1 ** 2 * (ZETA32 - 1) * exp(pi * Q1 ** 2 / 64) * ey
    R3 = 2 * Q1 ** 2 * ZETA32 * T * ey
    return (R1 + R2 + R3) / (sqrt(2) * exp(pi * sqrt(nu / 2)) * (1 - exp(-2 * pi * sqrt(nu / 2))))

def M(nu): return 2 * sqrt(2) * sinh(pi * sqrt(mpf(nu) / 2))

if __name__ == "__main__":
    KX = 3000
    b = beta_list(KX)
    B = [(-1) ** k * b[k] for k in range(KX + 1)]
    bh = [0] * (KX + 1); s = 0
    for k in range(KX + 1):
        s += B[k]; bh[k] = s
    print("beta_1..12:", b[1:13])
    print("T_inf = %s, T_R = %s" % (mp.nstr(T_inf, 8), mp.nstr(T_R, 8)))
    # agreement with exact values (consistency of the analysis)
    for k in [50, 100, 200, 500, 1000, 3000]:
        nu = k - mpf(1) / 8
        act = abs(bh[k] / (-M(nu)) - 1)
        print("  k=%d: |bhat/main - 1| = %s <= E = %s" % (k, mp.nstr(act, 4), mp.nstr(Ebeta(nu), 4)))
        assert act <= Ebeta(nu)
    # B_k < 0 for k >= K0 from the envelope: M(nu)(1-E(nu)) > M(nu-1)(1+E(nu-1))
    K0 = None
    for k in range(2, KX + 1):
        nu = k - mpf(1) / 8
        ok = M(nu) * (1 - Eplus(nu)) > M(nu - 1) * (1 + Eplus(nu - 1))
        if ok and K0 is None: K0 = k
        if not ok: K0 = None
    print("analytic: F_k(1) < 0 for all k >= K0 =", K0, "(checked up to", KX, "and the majorant decreases)")
    # beyond KX: log M(nu)/M(nu-1) >= pi/(2 sqrt(2 nu)) * 0.99 while log((1+E)/(1-E)) <= 3E, E sqrt(nu) decreasing
    for nu in [KX, 10 ** 4, 10 ** 5, 10 ** 6]:
        nu = mpf(nu)
        assert 3 * Eplus(nu - 1) < mpf('0.99') * pi / (2 * sqrt(2 * nu))
        assert Eplus(10 * nu) * sqrt(10 * nu) < Eplus(nu) * sqrt(nu)
    # exact range
    assert all(B[k] < 0 for k in range(7, KX + 1)) and all(B[k] < 0 for k in range(1, 8, 2))
    print("exact: F_k(1) < 0 for 7 <= k <= %d and k in {1,3,5}; F_2, F_4, F_6 =" % KX, B[2], B[4], B[6])
    print("so beta_k > 0 for every odd k and beta_k < 0 for every even k >= 8; beta_2, beta_4, beta_6 =", b[2], b[4], b[6])
    # growth: log|beta_k| - pi sqrt(k/2)
    for k in [100, 1000, 3000]:
        nu = k - mpf(1) / 8
        print("  k=%d: |beta_k| / (M(nu)-M(nu-1)) = %s" % (k, mp.nstr(abs(b[k]) / (M(nu) - M(nu - 1)), 12)))
    # explicit error for |beta_k| = (M(nu)-M(nu-1))(1 + theta eps'_k), eps'_k = (E(nu)M(nu)+E(nu-1)M(nu-1))/(M(nu)-M(nu-1))
    def epsp(k):
        nu = k - mpf(1) / 8
        return (Eplus(nu) * M(nu) + Eplus(nu - 1) * M(nu - 1)) / (M(nu) - M(nu - 1))
    ev = [(k, epsp(k)) for k in [80, 100, 200, 400, 1000, 3000, 10 ** 4]]
    print("eps'_k:", [(k, mp.nstr(v, 4)) for k, v in ev])
    assert all(ev[i + 1][1] < ev[i][1] for i in range(len(ev) - 1)) and ev[0][1] < 1
    print("BETA OK")
