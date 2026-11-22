"""Polar linear algebra for extremal N=2 elliptic genera.
For index m, the basis of weak Jacobi forms of weight 0 is
E4^a E6^b phi_{-2,1}^c phi_{0,1}^{m-c}, 4a+6b=2c.  We compute the polar
coefficients (l,n), 1<=l<=m, 4mn-l^2<0, and test whether p_ext is in the span."""
import numpy as np, sys
from jacobi import phi_m21, phi_01, E4, E6, ext_polar_source

def to_arr(F, X, L, p):
    """dict series -> array [X+1][2L+1] mod p (x exponents 0..X, l in -L..L)"""
    A = np.zeros((X+1, 2*L+1), dtype=np.int64)
    for k, v in F.items():
        if k <= X:
            for l, c in v.items():
                if -L <= l <= L: A[k, l+L] = c % p
    return A

def amul(A, B, p):
    X = A.shape[0]-1; L = (A.shape[1]-1)//2
    C = np.zeros_like(A)
    for i in range(X+1):
        if not A[i].any(): continue
        for j in range(X+1-i):
            if not B[j].any(): continue
            c = np.convolve(A[i], B[j]) % p     # length 4L+1, centre at 2L
            C[i+j] = (C[i+j] + c[L:3*L+1]) % p
    return C

def polar_points(m):
    pts = []
    for l in range(1, m+1):
        n = 0
        while 4*m*n - l*l < 0:
            pts.append((l, n)); n += 1
    return pts

def modforms(w):
    return [(a, (w-4*a)//6) for a in range(w//4+1) if (w-4*a) % 6 == 0]

def system(m, p):
    X = m//2 + 1           # x = q^(1/2); need n <= m/4
    L = m + X + 2          # truncation in l is safe: |l| > L never feeds |l| <= m? (checked below)
    # products only increase |l| by bounded amounts per factor, but large |l| from one factor can be
    # cancelled by the other; keep full range: index-m forms have l^2 <= 4mn+m^2 per term,
    # each factor of index 1 has l^2 <= 4n+1, so |l| <= m + 2n overall, n <= X/2.
    L = m + X + 2
    P21 = to_arr(phi_m21(X), X, L, p); P01 = to_arr(phi_01(X), X, L, p)
    e4 = to_arr(E4(X), X, L, p); e6 = to_arr(E6(X), X, L, p)
    pw21 = [None]*(m+1); pw01 = [None]*(m+1)
    I = np.zeros((X+1, 2*L+1), dtype=np.int64); I[0, L] = 1
    pw21[0] = I; pw01[0] = I
    for c in range(1, m+1):
        pw21[c] = amul(pw21[c-1], P21, p); pw01[c] = amul(pw01[c-1], P01, p)
    pe4 = [I]; pe6 = [I]
    for c in range(1, m//2+2):
        pe4.append(amul(pe4[-1], e4, p)); pe6.append(amul(pe6[-1], e6, p))
    pts = polar_points(m)
    cols = []
    for c in range(m+1):
        base = amul(pw21[c], pw01[m-c], p)
        for a, b in modforms(2*c):
            M = amul(amul(base, pe4[a], p), pe6[b], p)
            cols.append([M[2*n, l+L] for (l, n) in pts])
    A = np.array(cols, dtype=np.int64).T % p
    S = to_arr(ext_polar_source(m, X), X, L, p)
    bvec = np.array([S[2*n, l+L] for (l, n) in pts], dtype=np.int64) % p
    return A, bvec, pts

def rank_mod(M, p):
    M = M.copy() % p; r = 0; rows, cols = M.shape
    for c in range(cols):
        piv = None
        for i in range(r, rows):
            if M[i, c]: piv = i; break
        if piv is None: continue
        M[[r, piv]] = M[[piv, r]]
        inv = pow(int(M[r, c]), p-2, p)
        M[r] = (M[r]*inv) % p
        for i in range(rows):
            if i != r and M[i, c]:
                M[i] = (M[i] - M[i, c]*M[r]) % p
        r += 1
        if r == rows: break
    return r

if __name__ == "__main__":
    p = 33554393  # prime < 2^25
    lo, hi = int(sys.argv[1]), int(sys.argv[2])
    for m in range(lo, hi+1):
        A, b, pts = system(m, p)
        j = A.shape[1]; P = A.shape[0]
        rA = rank_mod(A, p); rAb = rank_mod(np.column_stack([A, b]), p)
        status = "consistent" if rAb == rA else "INCONSISTENT"
        print(f"m={m} P(m)={P} j(m)={j} rank(A)={rA} rank(A|b)={rAb} {status}", flush=True)
