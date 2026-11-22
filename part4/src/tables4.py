"""Data for the tables and figures of Part IV."""
from constraints import odd_S, lattice_form
from sseries import series, inv_series
from flint import fmpz_poly
from exist import solve
from math import isqrt, log, pi, sqrt

# odd levels
M = 4001
S = odd_S(M)
print("ODD TABLE (m, (-1)^((m-1)/2) S(m), C(2m-1), r(2m-1), 3m-1)")
for m in range(1, 62, 2):
    C, r = lattice_form(m); s = (-1)**((m-1)//2)
    print(m, s*S[m], C, r, 3*m-1)
sols = [m for m in range(1, M+1, 2) if (-1)**((m-1)//2)*S[m] == 3*m-1]
print("odd solutions m<=%d:" % M, sols)
# effective threshold
m = 1
while True:
    N = 2*m-1
    if not (pi*N/2 - 6*sqrt(N) - 4 > 3*m-1): last = m
    m += 2
    if m > 5000: break
print("largest odd m violating the analytic lower bound:", last)
print("ODDPLOT", " ".join("(%d,%.4f)" % (m, (-1)**((m-1)//2)*S[m]/m) for m in range(1, 801, 2)))
# even levels
Ssz = 200
G, T, D = series(Ssz)
def rat(num, k):
    return num.mul_low(inv_series(fmpz_poly([1, 1]).pow_trunc(k, Ssz+1), Ssz+1), Ssz+1)
r1 = rat(fmpz_poly([1, -1]), 1); r3 = rat(fmpz_poly([0, 2, -2]), 3)
LL = r1.mul_low(G.mul_low(T, Ssz+1), Ssz+1); R2 = r1.mul_low(D, Ssz+1); R1 = r3.mul_low(G, Ssz+1)
print("TWO TABLE (m, L(m), -12m)")
for m in range(2, 63, 4): print(m, int(LL[m//2]), -12*m)
print("ZERO TABLE (m, R1, R2, X)")
for m in range(4, 61, 4):
    N = m//2; X = -(4*m-2)*int(R2[N]) - int(R1[N]) - 2*m*m
    print(m, int(R1[N]), int(R2[N]), X)
print("TWOPLOT", " ".join("(%d,%.4f)" % (m, log(-int(LL[m//2]))/log(10)) for m in range(6, 400, 4)))
print("GSERIES first", [int(G[i]) for i in range(12)])
print("T first", [int(T[i]) for i in range(8)])
# extremal genera: q^1 row
for m in [1, 2, 3, 4, 5, 7, 8, 11, 13]:
    (sol, labels, phi), rank, J = solve(m)
    row1 = {l: phi.get(2, {}).get(l, 0) for l in range(-m-2, m+3) if phi.get(2, {}).get(l, 0) != 0}
    print("GENUS m=%d j=%d c(0,0)=%s c(1,l>=0)=%s" % (m, J, phi[0].get(0), [str(row1.get(l, 0)) for l in range(0, m+3)]))
