"""Key inequalities (L1) 2048 e_1(k-1) > |beta_k| and (L2) 2048 (e_1(k)-e_1(k-1)) > |beta_k|+|beta_{k+1}|,
exactly for 1 <= k <= KMAX, plus the analytic tail beyond KMAX."""
import time
from math import log, sqrt, pi, isqrt
from flint import fmpz_poly
from ramond import E_series
from superext import beta_list, dlist
KMAX = 3000
t0 = time.time()
beta = beta_list(KMAX + 2)
e1 = [int(v) for v in E_series(1, KMAX + 1).coeffs()]
print("e_1(0..8) =", e1[:9])
worst1 = worst2 = None
for k in range(1, KMAX + 1):
    lhs, rhs = 2048 * e1[k - 1], abs(beta[k])
    assert lhs > rhs, ("L1", k)
    lhs2, rhs2 = 2048 * (e1[k] - e1[k - 1]), abs(beta[k]) + abs(beta[k + 1])
    assert lhs2 > rhs2, ("L2", k)
    # record smallest ratio (in log10)
    r1 = log(lhs) - log(max(rhs, 1)); r2 = log(lhs2) - log(max(rhs2, 1))
    if worst1 is None or r1 < worst1[0]: worst1 = (r1, k)
    if worst2 is None or r2 < worst2[0]: worst2 = (r2, k)
print("L1, L2 hold for 1<=k<=%d" % KMAX)
print("min log ratio L1 = %.4f at k=%d ; L2 = %.4f at k=%d" % (worst1[0], worst1[1], worst2[0], worst2[1]))
for k in (1, 2, 3, 4, 8, 10, 100, 1000, 3000):
    print("k=%d: 2048 e_1(k-1)=%.4e  |beta_k|=%d..." % (k, 2048 * e1[k - 1], abs(beta[k])) if k < 20 else
          "k=%d: log10(2048 e_1(k-1))=%.2f  log10|beta_k|=%.2f" % (k, log(2048 * e1[k - 1], 10), log(abs(beta[k]), 10)))
# majorant |beta_k| <= e^{pi sqrt(k/2)} (1+24k^2) checked for all k <= KMAX
for k in range(1, KMAX + 1):
    assert log(abs(beta[k]) + 1e-300) <= pi * sqrt(k / 2) + log(1 + 24 * k * k), k
print("|beta_k| <= e^{pi sqrt(k/2)}(1+24k^2) checked for k<=%d" % KMAX)
# q_d lower bound and e_1 >= q_d(floor(N/23))^23, e_1(N)-e_1(N-1) >= q_d(floor(N/22))^22
qd = [0] * (KMAX + 2); qd[0] = 1
for part in range(1, KMAX + 2):
    for i in range(KMAX + 1, part - 1, -1): qd[i] += qd[i - part]
def Lstar(M):
    L = 0
    while (L + 1) * (L + 4) // 2 < M: L += 1
    return L
for M in range(0, KMAX + 2):
    assert qd[M] >= 2 ** Lstar(M), M
    assert Lstar(M) >= sqrt(2 * M + 2.25) - 2.5 - 1e-12, M
for N in range(0, KMAX + 1):
    assert e1[N] >= qd[N // 23] ** 23, N
    if N >= 1: assert e1[N] - e1[N - 1] >= qd[N // 22] ** 22, N
print("q_d(M) >= 2^{L*(M)}, L*(M) >= sqrt(2M+9/4)-5/2, e_1 lower bounds: checked to %d" % KMAX)
# analytic tail: phi1(k), phi2(k) > 0 for k > KMAX and increasing
def phi1(k):
    L = sqrt(2 * (k - 23) / 23) - 2.5
    return log(2048) + 23 * log(2) * L - pi * sqrt(k / 2) - log(1 + 24 * k * k)
def phi2(k):
    L = sqrt(2 * (k - 22) / 22) - 2.5
    return log(2048) + 22 * log(2) * L - log(2) - pi * sqrt((k + 1) / 2) - log(1 + 24 * (k + 1) ** 2)
print("phi1(3001)=%.3f  phi2(3001)=%.3f" % (phi1(KMAX + 1), phi2(KMAX + 1)))
for k in range(KMAX + 1, 200000, 7):
    assert phi1(k + 7) > phi1(k) > 0 and phi2(k + 7) > phi2(k) > 0
print("phi1, phi2 positive and increasing on sampled grid to 2e5 (derivative argument covers all k)")
print("time %.1fs" % (time.time() - t0))
