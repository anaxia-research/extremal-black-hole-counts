"""m = 0 mod 4 constraint derived directly.  Variable s = q^(1/2)."""
from fractions import Fraction
from collections import defaultdict
from flint import fmpz_poly

def sig(n): return sum(d for d in range(1, n+1) if n % d == 0)

def V2(m):
    """W(s,y) = (1-s^2) prod_{n>=1} (1+y s^{2n+1})(1+s^{2n+1}/y)/(1-s^{2n})^2 up to s^{m/2}.
    Returns dict s-exp -> dict y-exp -> int."""
    S = m//2
    W = {0: {0: 1}}
    def mul(A, B):
        C = defaultdict(lambda: defaultdict(int))
        for a, pa in A.items():
            for b, pb in B.items():
                if a+b > S: continue
                for l1, c1 in pa.items():
                    for l2, c2 in pb.items(): C[a+b][l1+l2] += c1*c2
        return {k: dict(v) for k, v in C.items()}
    W = mul(W, {0: {0: 1}, 2: {0: -1}})
    for n in range(1, S+1):
        e = 2*n+1
        if e <= S:
            W = mul(W, {0: {0: 1}, e: {1: 1}}); W = mul(W, {0: {0: 1}, e: {-1: 1}})
        if 2*n <= S:
            geo = {2*n*k: {0: k+1} for k in range(S//(2*n)+1)}   # 1/(1-s^{2n})^2
            W = mul(W, geo)
    return W

def theta_T(S):
    """T = theta4^4 - theta2^4 in s = q^(1/2), integer exponents (theta2^4 = 16 q^{1/2}(...))"""
    # theta4 = sum (-1)^n s^{n^2}; theta2^4: theta2 = sum q^{(n+1/2)^2/2} = sum s^{(n+1/2)^2}; theta2^4 exps = sum of 4 (n+1/2)^2 = integer
    X = 4*S + 4
    t4 = [0]*(X+1); n = 0
    while n*n <= S:
        t4[n*n] += 1 if n == 0 else 2*(-1)**n; n += 1
    t4 = fmpz_poly(t4[:S+1])
    # theta2 in u = s^(1/4): exponents (2n+1)^2
    t2 = [0]*(4*S+2); n = 0
    while (2*n+1)**2 <= 4*S+1:
        t2[(2*n+1)**2] += 2; n += 1
    t2 = fmpz_poly(t2)
    t24 = t2**4
    T = [int((t4**4)[i]) - int(t24[4*i]) for i in range(S+1)]
    return T

def X_constraint(m):
    S = m//2               # s-exponent of q^{m/4}
    W = V2(m)
    T = theta_T(S)
    # V(tau,0) = s^{-S} W(s,1): h(-a) for a = (S-j)/2 in s-units: coefficient at s^j, exponent j-S
    h = {}
    for j, pj in W.items():
        h[j - S] = sum(pj.values())
    # C0 = 12m - sum_{a>0} h(-a) T_a   (a in s-units: exponent -k, k>0)
    C0 = 12*m - sum(h.get(-k, 0)*T[k] for k in range(1, S+1))
    # E2 = 1 - 24 sum sigma(n) q^n = s-exps 2n
    E2part = sum(h.get(-2*n, 0)*(-24*sig(n)) for n in range(1, S//2+1))
    k2 = sum(l*l*c for l, c in W.get(S, {}).items())
    X = Fraction(k2) - Fraction(m, 6)*(E2part + C0)
    return X, C0, k2

if __name__ == "__main__":
    import sys
    M = int(sys.argv[1]) if len(sys.argv) > 1 else 80
    for m in range(4, M+1, 4):
        X, C0, k2 = X_constraint(m)
        print(m, "X =", X, " C0 =", C0)
