"""Exact s-series (s = q^{1/2}) for the constraints with m even, via flint."""
from flint import fmpz_poly
import time, sys
from math import log

def inv_series(f, n):
    """1/f mod s^n, f(0) = +-1, by Newton iteration"""
    c0 = int(f[0]); assert c0 in (1, -1)
    g = fmpz_poly([c0]); k = 1
    while k < n:
        k = min(2*k, n)
        e = f.mul_low(g, k)
        g = g.mul_low(fmpz_poly([2]) - e, k)
    return g

def P3inv(S):
    # prod (1-s^{2n})^3 = sum_k (-1)^k (2k+1) s^{k(k+1)}  (Jacobi)
    c = [0]*(S+1); k = 0
    while k*(k+1) <= S:
        c[k*(k+1)] += (-1)**k*(2*k+1); k += 1
    return inv_series(fmpz_poly(c), S+1)

def series(S):
    c3 = [0]*(S+1); c4 = [0]*(S+1); cd = [0]*(S+1); n = 0
    while n*n <= S:
        c3[n*n] += 1 if n == 0 else 2
        c4[n*n] += 1 if n == 0 else 2*(-1)**n
        cd[n*n] += n*n
        n += 1
    psi = [0]*(S+1); n = 0
    while n*(n+1) <= S:
        psi[n*(n+1)] += 2; n += 1
    t3, t4, dt3, ps = fmpz_poly(c3), fmpz_poly(c4), fmpz_poly(cd), fmpz_poly(psi)
    P = P3inv(S)
    G = t3.mul_low(P, S+1)
    t24 = fmpz_poly([0, 1]).mul_low(ps.pow_trunc(4, S+1), S+1)
    T = t4.pow_trunc(4, S+1) - t24
    D = dt3.mul_low(P, S+1)
    return G, T, D

def rat(S, num, den_pow):
    """series of num(s)/(1+s)^den_pow"""
    inv = inv_series(fmpz_poly([1, 1]).pow_trunc(den_pow, S+1), S+1)
    return num.mul_low(inv, S+1)

if __name__ == "__main__":
    S = int(sys.argv[1])
    t0 = time.time()
    G, T, D = series(S)
    r1 = rat(S, fmpz_poly([1, -1]), 1)            # (1-s)/(1+s)
    r3 = rat(S, fmpz_poly([0, 2, -2]), 3)        # 2s(1-s)/(1+s)^3
    LL = r1.mul_low(G.mul_low(T, S+1), S+1)
    R2 = r1.mul_low(D, S+1)
    R1 = r3.mul_low(G, S+1)
    print("built in %.1fs" % (time.time()-t0))
    bad2 = [2*N for N in range(1, S+1, 2) if not (int(LL[N]) < -24*N)]
    print("m=2 mod 4, m<=%d, failures of L < -12m:" % (2*S), bad2)
    bad0 = []
    for N in range(2, S+1, 2):
        m = 2*N
        X = -(4*m-2)*int(R2[N]) - int(R1[N]) - 2*m*m
        if not X < 0: bad0.append((m, X))
    print("m=0 mod 4, m<=%d, failures of X < 0:" % (2*S), bad0)
    negR2 = [2*N for N in range(2, S+1, 2) if int(R2[N]) <= 0]
    print("m=0 mod 4 with R2 <= 0:", negR2[:20])
    for N in [S//4, S//2, S-2]:
        N -= N % 2
        g = int(G[N]); 
        print(N, "log G_N=%.3f pi*sqrt(N)=%.3f  R2/G=%.5f  R1/G=%.5f  L/G=%.5f" % (log(g), 3.14159265*N**0.5, int(R2[N])/g, int(R1[N])/g, int(LL[N+1])/int(G[N+1])))
