"""Exact finite checks: cross-validation, positivity, monotonicity in n and k, log-concavity."""
import sys, time
from jcoeffs import jcoeffs
from extremal import p2list, pentagonal, primaries, direct_Z, hecke_w

L = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
K = int(sys.argv[2]) if len(sys.argv) > 2 else 40
t0 = time.time()
c = jcoeffs(L)
print(f"j coefficients to N={L}: c(1)={c[1]}, c(2)={c[2]}, c(3)={c[3]}")
p2 = p2list(K + L + 5)
pent = pentagonal(L + K + 5)

# 1. cross-check of the Hecke formula against direct polar-part matching
for k in range(1, 13):
    T = min(60, L // k)
    w = hecke_w(k, T, c, p2)
    Z = direct_Z(k, T, c)
    assert all(w[t] == Z[t] for t in range(-k, T + 1)), k
print("Hecke formula agrees with direct J-polynomial computation for k<=12, n<=60")

PR = {}
for k in range(1, K + 1):
    PR[k] = primaries(k, L // k, c, p2, pent)[0]
assert PR[1][1] == 196883 and PR[1][2] == 21296876 and PR[4][1] == 81026609426
print("Witten's values reproduced: 196883, 21296876 (k=1); 81026609426 (k=4, h=5)")

tot = sum(len(PR[k]) for k in PR)
neg = [(k, n) for k in PR for n in PR[k] if PR[k][n] <= 0]
print(f"positivity: {tot} values p_k(n), k<={K}, k*n<={L}; non-positive: {neg}")
mn = [(k, n) for k in PR for n in PR[k] if n + 1 in PR[k] and not PR[k][n + 1] > PR[k][n]]
print(f"strict increase in n: failures {mn}")
mk = [(k, n) for k in range(1, K) for n in PR[k + 1] if n in PR[k] and not PR[k + 1][n] > PR[k][n]]
print(f"strict increase in k: failures {mk}")
lc = [(k, n) for k in PR for n in PR[k] if n > 1 and n + 1 in PR[k] and not PR[k][n] ** 2 > PR[k][n - 1] * PR[k][n + 1]]
print(f"strict log-concavity in n: failures {lc}")
# the finite range used in the proof: k<=9, n<=9k and k<=10 for the k-monotonicity
req = [(k, n) for k in range(1, 10) for n in range(1, 9 * k + 1)]
assert all(PR[k][n] > 0 for k, n in req)
print(f"finite range of Theorem A (k<=9, n<=9k): {len(req)} values, all positive")
print("table p_k(1), p_k(2):")
for k in range(1, 11):
    print(k, PR[k][1], PR[k][2])
print(f"time {time.time()-t0:.1f}s")
