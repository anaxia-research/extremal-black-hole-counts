"""Light regime of Part II: k >= 10, 1 <= j <= 9k.  Certifies the constants of Section 'Light black holes'."""
from consts2 import *
from math import ceil as iceil
mp.dps = 40
LAM = mpf(9) / 10
cB = 2 * pi - pi / sqrt(2)          # exponent constant in the small-index terms

C1 = (1 + delta_bar(9)) * (1 + eta_bar(9)) / ((1 - delta_bar(9)) * (1 - eta_bar(9)))
B0 = (1 + eta_bar(1)) * (1 + delta_bar(1)) / ((1 - eta_bar(20)) * (1 - delta_bar(20)))

def rho_main(k):
    k = mpf(k)
    return C1 * (9 * k / (9 * k - 1)) ** (mpf(3) / 4) * exp(-2 * pi * sqrt(LAM * k) * (sqrt(9 * k) - sqrt(9 * k - 1)))

def tail(L):        # sum_{l >= L} e^{-cB sqrt l} <= e^{-c sqrt L} + (2/c^2)(1 + c sqrt L) e^{-c sqrt L}
    L = mpf(L)
    return exp(-cB * sqrt(L)) * (1 + 2 / cB ** 2 * (1 + cB * sqrt(L)))

def beta(k):
    return B0 * 2 ** (mpf(3) / 4) * tail(k // 10 + 1)

def pi_k(k):        # q(k+1)/w_k(1) bound, decreasing in k
    k = mpf(k)
    return exp(pi * sqrt((k + 1) / 2)) / (k * (1 - eta_bar(10)) * (1 - delta_bar(10)) * Alow(k))

def lower(r, p): return 1 - r - r ** 3 - r ** 15 / (1 - r) - p
def upper(r, p): return 1 + r ** 6 + r ** 10 + r ** 28 / (1 - r) + p

def C2(k):
    n = iceil(k / 2)
    return (1 - eta_bar(n)) * (1 - delta_bar(n)) / ((1 + eta_bar(n)) * (1 + delta_bar(n)))

def Cm(m):
    return (1 - eta_bar(m)) * (1 - delta_bar(m)) / ((1 + eta_bar(m)) * (1 + delta_bar(m)))

if __name__ == "__main__":
    print("C1 = %s, B0 = %s" % (mp.nstr(C1, 10), mp.nstr(B0, 8)))
    # (1) monotonicity in k of rho_main, beta, pi_k (formulas) -- grid confirmation
    ks = [10, 11, 20, 50, 100, 1000, 10 ** 4, 10 ** 6]
    assert all(rho_main(ks[i + 1]) < rho_main(ks[i]) for i in range(len(ks) - 1))
    assert all(beta(ks[i + 1]) <= beta(ks[i]) for i in range(len(ks) - 1))
    assert all(pi_k(ks[i + 1]) < pi_k(ks[i]) for i in range(len(ks) - 1))
    # t-monotonicity of (3/4) log(t/(t-1)) - 2 pi sqrt(m)(sqrt t - sqrt(t-1)) for m >= 9, t >= 2: pi sqrt(m)(t-1) > 1.5 sqrt(t)
    assert 3 * pi * 1 > mpf('1.5') * sqrt(2)
    rho = rho_main(10) + beta(10)
    print("rho_main(10) = %s, beta_10 = %s, rho = %s (limit %s)" % (mp.nstr(rho_main(10), 8), mp.nstr(beta(10), 6),
          mp.nstr(rho, 8), mp.nstr(C1 * exp(-pi * sqrt(LAM) / 3), 8)))
    p10 = pi_k(10)
    lo, up = lower(rho, p10), upper(rho, p10)
    print("pi_10 = %s; lower factor %s, upper factor %s, lower/rho = %s" % (mp.nstr(p10, 4), mp.nstr(lo, 8), mp.nstr(up, 8), mp.nstr(lo / rho, 6)))
    assert lo > 0 and lo / rho > up
    # (2) term-by-term ratio for the k-monotonicity lemma: C(m) e^{(1-2^{-1/2})/sqrt(m+1) * pi} >= 1 for all m
    assert all(Cm(m) * exp(pi * (1 - 1 / sqrt(2)) / sqrt(m + 1)) >= 1 for m in range(1, 3001))
    assert 4 * (delta_plus(3001) + eta_plus(3001)) < pi * (1 - 1 / sqrt(2)) / sqrt(10 ** 12)  # and both decay faster than m^{-1/2}
    print("C(m) e^{pi(1-2^-1/2)/sqrt(m+1)} >= 1 for m <= 3000 (min %s); C2(10) = %s" %
          (mp.nstr(min(Cm(m) * exp(pi * (1 - 1 / sqrt(2)) / sqrt(m + 1)) for m in range(1, 3001)), 6), mp.nstr(C2(10), 8)))
    # (3) k-monotonicity in the light range
    def r_alpha(alpha, k):
        return C1 * ((mpf(4) / 3) ** (mpf(3) / 4) * exp(-pi * sqrt(LAM / alpha)) + 2 ** (mpf(3) / 4) * exp(-pi * sqrt(LAM * k / 3))) + beta(k)
    def margin(alpha, k, a_lo=None):
        # left side with exponential at alpha_lo, r at alpha (worst ends)
        a_lo = alpha if a_lo is None else a_lo
        r = r_alpha(alpha, k)
        L = C2(k) * exp(pi * sqrt(a_lo * k / (k + 1))) * lower(r, pi_k(k + 1))
        return L / upper(r, pi_k(k))
    for (k0, amin, label) in [(10, mpf(1) / 2000, "10<=k<=2000"), (2000, mpf(1) / 2000, "k>2000, alpha>=1/2000")]:
        # j = 1 handled separately: rho = 0
        n = 18000
        worst = None
        for i in range(n):
            a0 = amin + (9 - amin) * i / n; a1 = amin + (9 - amin) * (i + 1) / n
            mrg = margin(a1, k0, a0)
            worst = mrg if worst is None or mrg < worst else worst
        print(label, ": min ratio of the two sides =", mp.nstr(worst, 6))
        assert worst > 1
    # j = 1 for all k >= 10:  C2(k) e^{pi/sqrt(k+1)} (1 - pi_{k+1}) > 1 + pi_k ;  check k <= 10^7 grid + asymptotics
    for k in [10, 100, 2000, 10 ** 4, 10 ** 5, 10 ** 6, 10 ** 7, 10 ** 9]:
        assert C2(k) * exp(pi / sqrt(k + 1)) * (1 - pi_k(k + 1)) > 1 + pi_k(k)
    # (4) case k > 2000, alpha < 1/2000:  log(N_{k+1}/N_k) >= g - 5 r1 - D_k,  D_k = 2(1-C2) + 5 beta_k + 2 pi_{k+1} + pi_k
    f2000 = 5 * C1 * 2 ** (mpf(3) / 4) / pi * sqrt(mpf('1.0005') * 2000) * exp(-pi * sqrt(LAM * 2000))
    print("f(2000) =", mp.nstr(f2000, 4))
    for K in [2000, 3000, 5000, 10 ** 4, 10 ** 5, 10 ** 6]:
        Dk = 2 * (1 - C2(K)) + 5 * beta(K) + 2 * pi_k(K + 1) + pi_k(K)
        rhs = (1 - f2000) * pi / sqrt(K + 1)
        print("  K=%d: D_k = %s < %s" % (K, mp.nstr(Dk, 4), mp.nstr(rhs, 4)))
        assert Dk < rhs
    print("LIGHT REGIME OK")
