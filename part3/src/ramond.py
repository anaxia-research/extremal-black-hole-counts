"""Part III: Ramond sector of Witten's N=1 extremal partition functions.
H_k(x) = (-1)^k f_k(-24 - 4096 x prod(1+x^n)^24)  (Ramond partition function, Witten (3.29)),
2 R_k(n) = sum_m (-1)^m h_{n-m^2}  (Ramond superconformal primaries).
Explicit positive formula:  sum_k H_k(x) p^k = D(-p) Lam(-p) / (1 - 4096 p x B(x) / A(-p))."""
import sys, time
from flint import fmpz_poly
from superext import kappa_list, dlist, ns_w, beta_list
from jcoeffs import jcoeffs

def trunc(P, L): return fmpz_poly([int(v) for v in P.coeffs()[:L + 1]])
def series_prod(factors, L):
    P = fmpz_poly([1])
    for f in factors: P = trunc(P * f, L)
    return P
def inv(P, L):  # power series inverse, P(0)=1
    c = [int(v) for v in P.coeffs()] + [0] * (L + 1); r = [0] * (L + 1); r[0] = 1
    for n in range(1, L + 1): r[n] = -sum(c[j] * r[n - j] for j in range(1, n + 1))
    return fmpz_poly(r)
def powtr(P, e, L):
    R = fmpz_poly([1]); B = P
    while e:
        if e & 1: R = trunc(R * B, L)
        B = trunc(B * B, L); e >>= 1
    return R
def mono(n, L, sgn=1):  # 1 + sgn x^n
    return fmpz_poly([1] + [0] * (n - 1) + [sgn]) if n <= L else fmpz_poly([1])

def Bser(L):   # prod (1+x^n)^24
    return powtr(series_prod([mono(n, L) for n in range(1, L + 1)], L), 24, L)

def H_direct(k, L, kap, c, d):
    """coefficients h_0..h_L of H_k by evaluating the polynomial f_k at g_R (independent route)."""
    xg = fmpz_poly([1, 0] + [kap[n] for n in range(1, 3 * k + 5)])
    pw = [fmpz_poly([1])]
    for r in range(1, k + 1): pw.append(trunc(pw[-1] * xg, 2 * k + 2))
    w = ns_w(k, 0, kap, c, d)
    rem = {t: w[t] for t in range(-k, 1)}; f = [0] * (k + 1)
    for r in range(k, -1, -1):
        cr = [int(v) for v in pw[r].coeffs()] + [0] * (3 * k + 5)
        a = rem[-r]; f[r] = a
        for t in range(-r, 1): rem[t] -= a * cr[t + r]
    gR = fmpz_poly([-24]) - 4096 * fmpz_poly([0, 1]) * Bser(L)
    H = fmpz_poly([0]); P = fmpz_poly([1])
    for r in range(k + 1):
        H += f[r] * P; P = trunc(P * gR, L)
    return [(-1) ** k * int(v) for v in H.coeffs()[:L + 1]] + [0] * (L + 1 - len(H.coeffs()[:L + 1]))

def E_series(r, K):
    """E_r(p) = (1+p) Lam(-p) prod (1-p^{2n-1})^{-(24r-1)} prod (1-p^{2n})^{-1}, r >= 1; r = 0: D(-p) Lam(-p)."""
    lam = [0] * (K + 1); lam[0] = 1
    for m in range(1, K + 1, 2):
        for j in range(m, K + 1, m): lam[j] += 24 * m     # 24 sum (2n-1) p^{2n-1}/(1-p^{2n-1})
    Lam = fmpz_poly(lam)
    odd = series_prod([mono(m, K, -1) for m in range(1, K + 1, 2)], K)   # prod (1-p^{2n-1})
    even = series_prod([mono(m, K, -1) for m in range(2, K + 1, 2)], K)  # prod (1-p^{2n})
    if r == 0:
        return trunc(trunc(fmpz_poly([1, 1]) * Lam, K) * trunc(odd * inv(even, K), K), K)
    return trunc(trunc(fmpz_poly([1, 1]) * Lam, K) * trunc(powtr(inv(odd, K), 24 * r - 1, K) * inv(even, K), K), K)

def rho_series(r, L):
    """rho_r(x) = prod (1+x^n)^{24r-3} * theta^+(x), theta^+ = sum x^{T_i}."""
    th = [0] * (L + 1); i = 0
    while i * (i + 1) // 2 <= L: th[i * (i + 1) // 2] = 1; i += 1
    return trunc(powtr(series_prod([mono(n, L) for n in range(1, L + 1)], L), 24 * r - 3, L) * fmpz_poly(th), L)

def R_explicit(K, L, beta):
    """R_k(n) for 1<=k<=K, 1<=n<=L from the positive formula."""
    e = {r: [int(v) for v in E_series(r, K).coeffs()] + [0] * (K + 2) for r in range(1, min(K, L) + 1)}
    rho = {r: [int(v) for v in rho_series(r, L).coeffs()] + [0] * (L + 2) for r in range(1, min(K, L) + 1)}
    R = {}
    for k in range(1, K + 1):
        for n in range(1, L + 1):
            s = sum(4096 ** r * rho[r][n - r] * e[r][k - r] for r in range(1, min(k, n) + 1))
            m = int(round(n ** 0.5))
            if m * m == n: s += 2 * (-1) ** m * beta[k]
            assert s % 2 == 0
            R[(k, n)] = s // 2
    return R, e

if __name__ == "__main__":
    t0 = time.time()
    kap = kappa_list(400); c = jcoeffs(120); d = dlist(500)
    beta = beta_list(60)
    # 1. independent check: R from H_direct vs explicit formula
    K, L = 14, 25
    R, e = R_explicit(K, L, beta)
    th = [0] * (L + 1)
    for m in range(-6, 7):
        if m * m <= L: th[m * m] += (-1) ** abs(m)
    bad = []
    for k in range(1, K + 1):
        h = H_direct(k, L, kap, c, d)
        assert h[0] == beta[k]
        for n in range(1, L + 1):
            two = sum(th[j] * h[n - j] for j in range(0, n + 1))
            if two != 2 * R[(k, n)]: bad.append((k, n))
            if h[n] % 4096: bad.append(('mod', k, n))
    print("explicit formula vs direct (k<=%d, n<=%d): mismatches" % (K, L), bad[:5])
    print("R_1(1..6) =", [R[(1, n)] for n in range(1, 7)])
    print("R_2(1..4) =", [R[(2, n)] for n in range(1, 5)])
    print("R_9(1..3) =", [R[(9, n)] for n in range(1, 4)])
    print("time %.1fs" % (time.time() - t0))
