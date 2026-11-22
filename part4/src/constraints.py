"""The three scalar constraints of GGKMO (5.19), (5.44) computed exactly.
Variable x = q^(1/8)."""
from flint import fmpz_poly, fmpz
from math import isqrt

def theta3(X):
    c = [0]*(X+1)
    n = 0
    while 4*n*n <= X:
        c[4*n*n] += 1 if n == 0 else 2; n += 1
    return fmpz_poly(c)

def theta4(X):
    c = [0]*(X+1); n = 0
    while 4*n*n <= X:
        c[4*n*n] += 1 if n == 0 else 2*(-1)**n; n += 1
    return fmpz_poly(c)

def theta2(X):
    c = [0]*(X+1); n = 0
    while (2*n+1)**2 <= X:
        c[(2*n+1)**2] += 2; n += 1
    return fmpz_poly(c)

def inv_series(f, X):
    # f(0) = 1
    g = fmpz_poly([1]); 
    # Newton iteration
    k = 1
    while k <= X:
        k = min(2*k, X+1)
        g = (g*(2 - (f*g).__class__([((f*g)[i]) for i in range(k)]))) 
        g = fmpz_poly([g[i] for i in range(k)])
    return g

def trunc(f, X):
    return fmpz_poly([f[i] for i in range(min(f.degree()+1, X+1))])

def mul(a, b, X):
    return a.mullow(b, X+1) if hasattr(a, 'mullow') else trunc(a*b, X)

def eta3inv(X):
    """1/prod(1-x^{8n})^3 = q^{1/8}/eta^3"""
    # prod(1-x^{8n})^3 = sum (-1)^k (2k+1) x^{8 k(k+1)/2}  (Jacobi)
    c = [0]*(X+1); k = 0
    while 4*k*(k+1) <= X:
        c[4*k*(k+1)] += (-1)**k*(2*k+1); k += 1
    return inv_series(fmpz_poly(c), X)

def ratio(X, step=4):
    """(1-x^4)/(1+x^4) = 1 + 2 sum_{k>=1} (-1)^k x^{4k}"""
    c = [0]*(X+1); c[0] = 1
    for k in range(1, X//step+1): c[step*k] = 2*(-1)**k
    return fmpz_poly(c)

def odd_S(M):
    """S(m) = [x^{2m-1}] (1-x^4)/(1+x^4) theta2 theta4, for odd m <= M"""
    X = 2*M
    F = mul(mul(ratio(X), theta2(X), X), theta4(X), X)
    return {m: int(F[2*m-1]) for m in range(1, M+1, 2)}

def lattice_form(m):
    N = 2*m-1; C = 0; r = 0
    b = 0
    for b in range(-isqrt(N//4)-1, isqrt(N//4)+2):
        R = N - 4*b*b
        if R < 0: continue
        s = isqrt(R)
        # odd u with u^2 <= R: u in {-s..s} odd
        k = (s+1)//2   # number of positive odd u <= s
        C += 2*k
        if s*s == R and s % 2 == 1: r += 2*(-1)**b
    return C, r

def even_L(M):
    """[x^{2m}] (1-x^4)/(1+x^4) theta3 (theta4^4 - theta2^4)/prod(1-x^{8n})^3 for m = 2 mod 4"""
    X = 2*M
    t3, t4, t2 = theta3(X), theta4(X), theta2(X)
    t44 = mul(mul(t4, t4, X), mul(t4, t4, X), X)
    t24 = mul(mul(t2, t2, X), mul(t2, t2, X), X)
    F = mul(mul(ratio(X), t3, X), mul(t44 - t24, eta3inv(X), X), X)
    return {m: int(F[2*m]) for m in range(2, M+1, 4)}

def R12(M):
    """R1 = 2[q^{m/4-1/8}] q^{1/2}(1-q^{1/2})^4/(1-q)^3 theta3/eta^3,
       R2 = [q^{m/4-1/8}] (1-q^{1/2})^2/(1-q) q d/dq theta3/eta^3, for m = 0 mod 4.
       In x: q^{m/4-1/8}/eta^3 -> coefficient of x^{2m} in (...) * q^{1/8}/eta^3."""
    X = 2*M
    t3 = theta3(X)
    e = eta3inv(X)
    # q d/dq theta3 = sum n^2/2 q^{n^2/2}: coefficient of x^{4n^2} is n^2/2 (times 2 for +-n)
    c = [0]*(X+1); n = 1
    while 4*n*n <= X:
        c[4*n*n] = n*n; n += 1      # 2 * n^2/2
    dt3 = fmpz_poly(c)
    # (1-x^4)^4/(1-x^8)^3 and (1-x^4)^2/(1-x^8)
    one_m = fmpz_poly([1,0,0,0,-1])
    den = [0]*(X+1)
    for k in range(0, X//8+1): den[8*k] = 1
    den = fmpz_poly(den)    # 1/(1-x^8)
    A1 = mul(mul(mul(one_m**4, den, X), mul(den, den, X), X), mul(t3, e, X), X)
    A1 = fmpz_poly([0,0,0,0] + [A1[i] for i in range(X-3)])   # times q^{1/2} = x^4
    A2 = mul(mul(one_m**2, den, X), mul(dt3, e, X), X)
    return {m: (2*int(A1[2*m]), int(A2[2*m])) for m in range(4, M+1, 4)}

if __name__ == "__main__":
    M = 400
    S = odd_S(M)
    bad = 0
    for m in range(1, M+1, 2):
        C, r = lattice_form(m)
        s = (-1)**((m-1)//2)
        if s*S[m] != 2*C - s*r: bad += 1
    print("lattice identity mismatches (odd m<=%d):" % M, bad)
    print("odd m with (-1)^((m-1)/2) S(m) = 3m-1:", [m for m in range(1, M+1, 2) if (-1)**((m-1)//2)*S[m] == 3*m-1])
    E = even_L(M)
    print("m=2 mod 4, L vs -12m:", [(m, E[m]) for m in list(E)[:8]])
    print("m=2 mod 4 with L = -12m:", [m for m in E if E[m] == -12*m])
    R = R12(M)
    print("m=0 mod 4 (R1,R2):", [(m,)+R[m] for m in list(R)[:10]])
    print("m=0 mod 4 satisfying -2m^2-R1-(4m-2)R2 = 0:", [m for m in R if -2*m*m - R[m][0] - (4*m-2)*R[m][1] == 0])
