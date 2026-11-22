"""Check: for every m <= 60, the coefficients of chi_NS at all exponents h <= 0
(except the D = 0, l = -m term) predicted from p_ext by spectral flow and charge
conjugation coincide with those of the NS vacuum character V."""
from jacobi import ext_polar_source, mul, factor, qmul, qpow, inv_q, scal
from fractions import Fraction
from collections import defaultdict

def V_ns(m):
    """V = q^{-m/4} (1-q) prod (1+y q^{n+1/2})(1+y^{-1} q^{n+1/2})/(1-q^n)^2, keys: s-exponent of q^{h+m/4} (s=q^{1/2})"""
    X = m//2 + 1
    F = {0: {0: 1}, 2: {0: -1}}
    for n in range(1, X+1):
        e = 2*n+1
        if e <= X:
            F = mul(F, factor(e, 1, 1, X), X); F = mul(F, factor(e, -1, 1, X), X)
    d = {0: 1}
    for n in range(1, X//2+1): d = qmul(d, {0: 1, 2*n: -1}, X)
    F = mul(F, scal(inv_q(qpow(d, 2, X), X)), X)
    return F, X

def check(m):
    X = m//2 + 1
    S = ext_polar_source(m, X)
    def pext(l, n):     # p_ext coefficient on P(m)
        if not (1 <= l <= m and n >= 0 and 4*m*n - l*l < 0): return None
        return S.get(2*n, {}).get(l, 0)
    V, XV = V_ns(m)
    bad = 0; count = 0
    # enumerate (n,l) with exponent h = n + l/2 + m/4 <= 0, n >= 0
    for n in range(0, m+1):
        for l in range(-3*m, m+1):
            h2 = 4*n + 2*l + m          # 4h
            if h2 > 0: continue
            D = 4*m*n - l*l
            if D >= 0:
                assert D == 0 and l == -m, (m, n, l)
                continue
            r = ((l + m) % (2*m)) - m     # r in [-m, m)
            if r == -m: r = m
            n2 = Fraction(D + r*r, 4*m)
            if n2.denominator != 1 or n2 < 0 or D < -m*m: pred = 0
            else:
                pred = pext(abs(r), int(n2)); pred = 0 if pred is None else pred
            pred_ns = (-1)**((l+m) % 2) * pred
            # V coefficient: s-exponent = 2(h + m/4) = 2n + l + m, y-exponent k = l + m
            vs = V.get(2*n + l + m, {}).get(l + m, 0)
            count += 1
            if vs != pred_ns: bad += 1
    return bad, count

for m in range(1, 61):
    b, c = check(m)
    if b: print("m=%d mismatches %d of %d" % (m, b, c))
print("checked m=1..60")
