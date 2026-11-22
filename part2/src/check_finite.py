"""Exact checks for Part II (N=1 extremal partition functions)."""
import sys, time
sys.path.insert(0, 'src'); sys.path.insert(0, '.')
from flint import fmpz_poly
from superext import *

t0 = time.time()
LJ = int(sys.argv[1]) if len(sys.argv) > 1 else 30000
kap = kappa_list(LJ); c = jcoeffs(LJ // 4 + 5); d = dlist(LJ + 100); tri = triangular(LJ + 100)
print("kappa(1..8) =", kap[1:9])

# 1. Replication (Hecke) formula against direct Faber polynomials, m <= 24, n <= 40
def faber_direct(M, L):
    xg = fmpz_poly([1, 0] + [kap[n] for n in range(1, L + M + 3)])   # q*g
    pw = [fmpz_poly([1])]
    for r in range(1, M + 1):
        pw.append(fmpz_poly([int(v) for v in (pw[-1] * xg).coeffs()[:L + 2 * M + 3]]))
    polys = []
    for m in range(M + 1):
        cm = [int(v) for v in pw[m].coeffs()] + [0] * (L + 2 * M + 5)
        D = {j - m: cm[j] for j in range(0, L + m + 1)}
        for r in range(m - 1, -1, -1):
            a = D.get(-r, 0)
            if a:
                for n, v in polys[r].items():
                    D[n] = D.get(n, 0) - a * v
        D = {n: v for n, v in D.items() if n <= L}
        if m > 0: D[0] = 0
        polys.append(D)
    return polys
P = faber_direct(24, 40)
bad = [(m, n) for m in range(1, 25) for n in range(1, 41) if P[m][n] != kappa_m(m, n, kap, c)]
print("replication formula vs direct Faber polynomials (m<=24, n<=40): mismatches", len(bad), bad[:5])
assert not bad

# 2. Witten's values
N9, _ = ns_primaries(9, 2, kap, c, d, tri)
N4, _ = ns_primaries(4, 2, kap, c, d, tri)
N1, _ = ns_primaries(1, 2, kap, c, d, tri)
print("k*=9, j=1:", N9[1], "(Witten 135149371);  k*=4:", N4[1], N4[2], "; k*=1:", N1[1], N1[2])
assert N9[1] == 135149371

# 3. Witten's table of beta
b = beta_list(400)
print("beta_1..10 =", b[1:11])
assert b[1:11] == [24, 24, 95, 1, 143, 1, 262, -213, 453, -261]
# independent check of beta for k <= 12: beta_k = (-1)^k F_k(R cusp) computed from Phi_k as polynomial in g, g(R)=-24
def beta_direct(k):
    # F_k = sum_{m} d(k-m) g_m + d(k); g_m(R) = value of Faber polynomial at g=-24: use P_m(g) coefficients
    xg = fmpz_poly([1, 0] + [kap[n] for n in range(1, 3 * k + 5)])
    # build polynomial f with F = f(g) by matching polar part
    pw = [fmpz_poly([1])]
    for r in range(1, k + 1): pw.append(fmpz_poly([int(v) for v in (pw[-1] * xg).coeffs()[:2 * k + 3]]))
    w = ns_w(k, 0, kap, c, d)
    rem = {t: w[t] for t in range(-k, 1)}
    f = [0] * (k + 1)
    for r in range(k, -1, -1):
        cr = [int(v) for v in pw[r].coeffs()] + [0] * (3 * k + 5)
        a = rem[-r]; f[r] = a
        for t in range(-r, 1): rem[t] -= a * cr[t + r]
    return (-1) ** k * sum(f[r] * (-24) ** r for r in range(k + 1))
assert all(beta_direct(k) == b[k] for k in range(1, 13))
print("beta generating function agrees with direct evaluation for k <= 12")
print("odd k<=400 with beta<=0:", [k for k in range(1, 401, 2) if b[k] <= 0],
      "; even 8<=k<=400 with beta>=0:", [k for k in range(8, 401, 2) if b[k] >= 0])

# 4. positivity and monotonicity of N_k(j), k <= 40, k*j <= LJ
cnt = 0; neg = []; nmono = []; kmono = []; lc = []
prev = None
NS = {}
for k in range(1, 41):
    T = LJ // k
    N, _ = ns_primaries(k, T, kap, c, d, tri)
    NS[k] = N
    cnt += T
    neg += [(k, j) for j in N if N[j] <= 0]
    nmono += [(k, j) for j in range(1, T) if N[j + 1] <= N[j]]
    lc += [(k, j) for j in range(2, T) if N[j] ** 2 <= N[j - 1] * N[j + 1]]
for k in range(1, 40):
    T = min(LJ // (k + 1), LJ // k)
    kmono += [(k, j) for j in range(1, T + 1) if NS[k + 1][j] <= NS[k][j]]
print("values:", cnt, "; non-positive:", neg[:10], "; not increasing in j:", nmono[:10],
      "; not increasing in k:", kmono[:10], "; not log-concave:", lc[:10])
print("table N_k(1), N_k(2):")
for k in range(1, 13):
    print(k, NS[k][1], NS[k][2])
print("time %.0fs" % (time.time() - t0))
