"""Exact rational extremal elliptic genera for m in {1,2,3,4,5,7,8,11,13};
exact nonexistence certificates (left-kernel vector) for m in {12,19,31,41}."""
from flint import fmpq_mat, fmpq, fmpz_mat
from jacobi import phi_m21, phi_01, E4, E6, ext_polar_source, mul
from polar import polar_points, modforms
import sys

def series_pow(F, k, X):
    R = {0: {0: 1}}
    for _ in range(k): R = mul(R, F, X)
    return R

def exact_system(m):
    X = m//2 + 1
    P21, P01, e4, e6 = phi_m21(X), phi_01(X), E4(X), E6(X)
    pw21 = [{0: {0: 1}}]; pw01 = [{0: {0: 1}}]
    for c in range(1, m+1):
        pw21.append(mul(pw21[-1], P21, X)); pw01.append(mul(pw01[-1], P01, X))
    pe4 = [{0: {0: 1}}]; pe6 = [{0: {0: 1}}]
    for c in range(1, m//2+2):
        pe4.append(mul(pe4[-1], e4, X)); pe6.append(mul(pe6[-1], e6, X))
    pts = polar_points(m)
    cols, labels, forms = [], [], []
    for c in range(m+1):
        base = mul(pw21[c], pw01[m-c], X)
        for a, b in modforms(2*c):
            M = mul(mul(base, pe4[a], X), pe6[b], X)
            cols.append([M.get(2*n, {}).get(l, 0) for (l, n) in pts]); labels.append((a, b, c)); forms.append(M)
    S = ext_polar_source(m, X)
    bvec = [S.get(2*n, {}).get(l, 0) for (l, n) in pts]
    return cols, bvec, pts, labels, forms, X

def solve(m):
    cols, b, pts, labels, forms, X = exact_system(m)
    P, J = len(pts), len(cols)
    A = fmpq_mat(P, J, [cols[j][i] for i in range(P) for j in range(J)])
    B = fmpq_mat(P, 1, b)
    # least squares via normal equations is exact if consistent; use rref of [A|B]
    Ab = fmpq_mat(P, J+1, [ (cols[j][i] if j < J else b[i]) for i in range(P) for j in range(J+1)])
    R, rank = Ab.rref()
    # find solution
    sol = [fmpq(0)]*J
    piv_row = 0
    pivots = []
    for i in range(rank):
        for j in range(J+1):
            if R[i, j] != 0:
                pivots.append(j); break
    if J in pivots: return None, rank, J
    for i, j in enumerate(pivots): sol[j] = R[i, J]
    # check
    for i in range(P):
        assert sum(sol[j]*cols[j][i] for j in range(J)) == b[i]
    # assemble the elliptic genus (truncated)
    phi = {}
    for j in range(J):
        if sol[j] == 0: continue
        for k, v in forms[j].items():
            for l, c in v.items():
                phi.setdefault(k, {}); phi[k][l] = phi[k].get(l, 0) + sol[j]*c
    return (sol, labels, phi), rank, J

if __name__ == "__main__":
    for m in [1, 2, 3, 4, 5, 7, 8, 11, 13, 12, 19]:
        res, rank, J = solve(m)
        if res is None:
            print(f"m={m}: NO rational solution (rank[A|b] includes the b column), j(m)={J}")
            continue
        sol, labels, phi = res
        integral = all(c.q == 1 for v in phi.values() for c in v.values())
        c0 = {l: phi[0].get(l, 0) for l in range(-m, m+1) if phi[0].get(l, 0) != 0}
        print(f"m={m}: solution found, j(m)={J}, all computed coefficients integral: {integral}; q^0 part {c0}")
