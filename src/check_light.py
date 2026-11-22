"""Light regime (k >= 10, n <= 9k): constants of Section 4 and the comparison inequalities."""
from mpmath import mp, mpf, sqrt, exp, log, pi, ceil
from bounds import *
from jcoeffs import jcoeffs

print("S = zeta(3/2)^2 - 1 =", mp.nstr(S, 10), "; delta constant =", mp.nstr(DELTA_C, 10))
# 1. sanity: delta bound against exact coefficients
c = jcoeffs(400)
mp.dps = 120
worst = max(abs(c[N] / A(N) - 1) / delta(N) for N in range(1, 401))
mp.dps = 40
print("max |c(N)/A(N)-1| / bound over N<=400:", mp.nstr(worst, 6), "(must be <= 1)")
assert worst <= 1
print("delta(1) =", mp.nstr(delta(1), 6), " delta(5) =", mp.nstr(delta(5), 6))

# 2. ratio bound rho for k >= 10, 2 <= t <= 9k
C1 = (1 + delta(5)) / (1 - delta(5)) * (1 + eta_div(5))
big4 = C1 * (mpf(4) / 3) ** (mpf(3) / 4) * exp(-2 * pi * sqrt(mpf(1) / 18))
big23 = C1 * 2 ** (mpf(3) / 4) * exp(-2 * pi * sqrt(mpf(10) / 6))
sm = max(small_part(k, 1) for k in range(10, 2001))   # t-1 >= 1, w(t) with t >= 2
print("C1 =", mp.nstr(C1, 10))
print("big ratio t>=4:", mp.nstr(big4, 8), "; t in {2,3}:", mp.nstr(big23, 8))
print("small-m part (max over 10<=k<=2000, decreasing beyond):", mp.nstr(sm, 8))
rho = big4 + sm
print("rho <=", mp.nstr(rho, 8))
negsum = rho + rho ** 2 + rho ** 12 / (1 - rho)
possum = rho ** 5 + rho ** 7 + rho ** 22 / (1 - rho)
print("lower factor 1 - (rho + rho^2 + rho^12/(1-rho)) =", mp.nstr(1 - negsum, 8))
print("upper factor 1 + rho^5 + rho^7 + rho^22/(1-rho) =", mp.nstr(1 + possum, 8))
# polar correction p2(k+n) / (k c(k n)) <= p(k+n)/(k c(kn)): bounded via p(k+n) <= p(10k) <= exp(pi sqrt(20k/3))
pol = max(p_bound(10 * k) / (k * A(k) * (1 - delta(k))) for k in range(10, 2001))
print("polar correction <=", mp.nstr(pol, 6))
low = 1 - negsum - pol
print("=> p_k(n) >= ", mp.nstr(low, 6), "* w_k(n) for k>=10, 1<=n<=9k")
# monotone in n
print("monotone in n: need low/rho > upper:", mp.nstr(low / rho, 6), ">", mp.nstr(1 + possum, 6))
assert low / rho > 1 + possum

# 3. decreasing tail of the small-m bound for k >= 2000 (log-derivative of the closed form)
#    f(k) = pi sqrt(2k/3) + log(k/2 (1+log(k/2))) - log k - 4 pi (1-1/sqrt2) sqrt(2k)  -> f'(k) < 0
import mpmath
fp = lambda k: pi / sqrt(6 * k) + 1 / (k * (1 + log(k / 2))) - 4 * pi * (1 - 1 / sqrt(2)) * sqrt(2) / (2 * sqrt(k))
print("derivative of log small-m bound at k=10:", mp.nstr(fp(mpf(10)), 6), "(negative and increasing to 0^-)")

# 4. monotonicity in k, light regime.
def rho_alpha(alpha, smv):
    return C1 * ((mpf(4) / 3) ** (mpf(3) / 4) * exp(-pi * sqrt(2 / alpha)) + 2 ** (mpf(3) / 4) * exp(-pi * sqrt(mpf(20) / 3))) + smv

def C2_of(k0):
    """Lower bound for the termwise ratio constant, valid for every k >= k0 (all pieces decrease in k)."""
    smk = max(small_part_same(k, 1) for k in range(k0, k0 + 200))   # decreasing: max at k0
    N0 = mpf(k0) / 2
    return (1 - smk) * (1 - delta(N0)) / ((1 + eta_div(N0)) * (1 + delta(N0)))

def alpha_grid_min(k0, smv, pol, C2, amin):
    worst = None
    n_steps = 18000
    grid = [amin + (9 - amin) * mpf(i) / n_steps for i in range(n_steps + 1)]
    for i in range(n_steps):
        a0, a1 = grid[i], grid[i + 1]
        r = rho_alpha(a1, smv)
        lhs = C2 * exp(2 * pi * sqrt(a0 * mpf(k0) / (k0 + 1))) * (1 - (r + r ** 2 + r ** 12 / (1 - r)) - pol)
        rhs = 1 + r ** 5 + r ** 7 + r ** 22 / (1 - r)
        g = lhs / rhs
        if worst is None or g < worst[0]:
            worst = (g, a0)
    return worst

# (a) 10 <= k <= 2000: alpha = n/k >= 1/2000
C2a = C2_of(10)
wa = alpha_grid_min(10, sm, pol, C2a, mpf(1) / 2000)
print("(a) 10<=k<=2000: C2 =", mp.nstr(C2a, 8), "; min LHS/RHS on alpha in [1/2000, 9]:", mp.nstr(wa[0], 8), "at", mp.nstr(wa[1], 5))
assert wa[0] > 1
# (b) k > 2000
C2b = C2_of(2000)
smb = small_part(2000, 1)
polb = p_bound(10 * 2000) / (2000 * A(2000) * (1 - delta(2000)))
wb = alpha_grid_min(2000, smb, polb, C2b, mpf(1) / 2000)
print("(b) k>2000, alpha>=1/2000: 1-C2 =", mp.nstr(1 - C2b, 6), "; min LHS/RHS:", mp.nstr(wb[0], 8))
assert wb[0] > 1
# (b') k > 2000, alpha = n/k < 1/2000.
#   log(p_{k+1}(n)/p_k(n)) >= g - D_k - 4 r_1,  g = 2 pi sqrt(n/(k+1)) >= 2 pi/sqrt(k+1),
#   r_1 = C1 2^{3/4} exp(-pi sqrt(2k/n)),  D_k = 2(1-C2(k)) + 4 beta_k + 2 pi_{k+1}.
#   4 r_1 / g <= f(x) := (4 C1 2^{3/4}/(2 pi)) sqrt(1.0005 x) exp(-pi sqrt(2x)),  x = k/n > 2000, f decreasing.
#   Need D_k < (1 - f(2000)) * 2 pi / sqrt(k+1); D_k sqrt(k+1) decreases (each piece is poly(k) exp(-1.1 sqrt k)).
with mp.workdps(400):
    f2000 = 4 * C1 * 2 ** (mpf(3) / 4) / (2 * pi) * sqrt(mpf(1.0005) * 2000) * exp(-pi * sqrt(mpf(4000)))
    print("(b') f(2000) = 4 r_1/g bound:", mp.nstr(f2000, 4))
    assert f2000 < mpf(10) ** -80
    for K in [2000, 3000, 5000, 10**4]:
        C2K = C2_of(K)
        DK = 2 * (1 - C2K) + 4 * small_part(K, 1) + 2 * p_bound(10 * (K + 1)) / ((K + 1) * A(K + 1) * (1 - delta(K + 1)))
        gain = (1 - f2000) * 2 * pi / sqrt(mpf(K + 1))
        print("   k =", K, ": D_k =", mp.nstr(DK, 4), " D_k sqrt(k+1) =", mp.nstr(DK * sqrt(mpf(K + 1)), 4), " < gain", mp.nstr(gain, 4))
        assert DK < gain
