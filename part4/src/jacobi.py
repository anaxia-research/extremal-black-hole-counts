"""Weak Jacobi forms of weight 0 as truncated double series.
A series is a dict: n (int or half-int via x=q^(1/2) exponent) -> dict l -> int.
We use the variable x = q^(1/2); keys are integer exponents of x."""
from collections import defaultdict
from fractions import Fraction

def mul(A, B, X):
    """multiply two series truncated at x-exponent <= X"""
    C = defaultdict(lambda: defaultdict(int))
    for a, pa in A.items():
        if a > X: continue
        for b, pb in B.items():
            if a + b > X: continue
            Cab = C[a+b]
            for l1, c1 in pa.items():
                for l2, c2 in pb.items():
                    Cab[l1+l2] += c1*c2
    return {k: {l: c for l, c in v.items() if c} for k, v in C.items()}

def one():
    return {0: {0: 1}}

def scal(s):  # q-only series from dict xexp->coef
    return {k: {0: c} for k, c in s.items()}

def inv_q(s, X):
    """inverse of a q-only power series (dict xexp->int) with s[0]=1"""
    r = {0: 1}
    for n in range(1, X+1):
        r[n] = -sum(s.get(k, 0)*r[n-k] for k in range(1, n+1))
    return r

def factor(e, l, sign, X):
    """1 + sign * x^e y^l"""
    return {0: {0: 1}, e: {l: sign}} if e <= X else one()

def phi_m21(X):
    """phi_{-2,1} = (y-2+1/y) prod (1-q^n y)^2 (1-q^n/y)^2/(1-q^n)^4, x=q^(1/2)"""
    F = {0: {1: 1, 0: -2, -1: 1}}
    for n in range(1, X//2+1):
        e = 2*n
        for _ in range(2):
            F = mul(F, factor(e, 1, -1, X), X); F = mul(F, factor(e, -1, -1, X), X)
    d = {0: 1}
    for n in range(1, X//2+1):
        d = qmul(d, {0: 1, 2*n: -1}, X)
    d4 = qpow(d, 4, X)
    F = mul(F, scal(inv_q(d4, X)), X)
    return F

def qmul(a, b, X):
    c = defaultdict(int)
    for i, u in a.items():
        for j, v in b.items():
            if i+j <= X: c[i+j] += u*v
    return dict(c)

def qpow(a, k, X):
    r = {0: 1}
    for _ in range(k): r = qmul(r, a, X)
    return r

def theta_ratio(kind, X):
    """theta_i(tau,z)^2/theta_i(tau)^2 for i=3 (kind=3), 4 (kind=4), times 4;
    for i=2 returns 4*theta_2(z)^2/theta_2^2 = (y+2+1/y) prod (1+q^n y)^2(1+q^n/y)^2/(1+q^n)^4"""
    if kind == 2:
        F = {0: {1: 1, 0: 2, -1: 1}}
        for n in range(1, X//2+1):
            for _ in range(2):
                F = mul(F, factor(2*n, 1, 1, X), X); F = mul(F, factor(2*n, -1, 1, X), X)
        d = {0: 1}
        for n in range(1, X//2+1): d = qmul(d, {0: 1, 2*n: 1}, X)
        F = mul(F, scal(inv_q(qpow(d, 4, X), X)), X)
        return F
    s = 1 if kind == 3 else -1
    F = {0: {0: 4}}
    for n in range(1, X+1):
        e = 2*n-1
        if e > X: break
        for _ in range(2):
            F = mul(F, factor(e, 1, s, X), X); F = mul(F, factor(e, -1, s, X), X)
    d = {0: 1}
    for n in range(1, X+1):
        if 2*n-1 > X: break
        d = qmul(d, {0: 1, 2*n-1: s}, X)
    F = mul(F, scal(inv_q(qpow(d, 4, X), X)), X)
    return F

def add(*As):
    C = defaultdict(lambda: defaultdict(int))
    for A in As:
        for k, v in A.items():
            for l, c in v.items(): C[k][l] += c
    return {k: {l: c for l, c in v.items() if c} for k, v in C.items() if any(v.values())}

def phi_01(X):
    return add(theta_ratio(2, X), theta_ratio(3, X), theta_ratio(4, X))

def sigma(n, k):
    return sum(d**k for d in range(1, n+1) if n % d == 0)

def E4(X):
    return scal({2*n: (240*sigma(n, 3) if n else 1) for n in range(X//2+1)})

def E6(X):
    return scal({2*n: (-504*sigma(n, 5) if n else 1) for n in range(X//2+1)})

def ext_polar_source(m, X):
    """(-1)^m SF_{1/2} chi_vac = (1-q) y^m prod (1-y q^{n+1})(1-y^{-1} q^n)/(1-q^n)^2"""
    F = {0: {m: 1}, 2: {m: -1}} if X >= 2 else {0: {m: 1}}
    for n in range(1, X//2+1):
        F = mul(F, factor(2*(n+1), 1, -1, X), X)
        F = mul(F, factor(2*n, -1, -1, X), X)
    d = {0: 1}
    for n in range(1, X//2+1): d = qmul(d, {0: 1, 2*n: -1}, X)
    F = mul(F, scal(inv_q(qpow(d, 2, X), X)), X)
    return F
