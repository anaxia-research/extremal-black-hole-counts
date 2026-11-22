from mpmath import mp, mpc, mpf, exp, pi
from exist import solve
mp.dps = 30
for m in [1, 2, 3, 4, 5, 7, 8, 11, 13]:
    (sol, labels, phi), rank, J = solve(m)
    nz = [(lab, s) for lab, s in zip(labels, sol) if s != 0]
    print("m=%2d  integral basis coefficients: %s  denominators: %s" % (m, all(s.q == 1 for _, s in nz), sorted(set(int(s.q) for _, s in nz))))
# kappa, K-hat at the Ramond cusp and theta identities
def th(i, t, K=80):
    s = exp(pi*1j*t)
    if i == 3: return 1 + 2*sum(s**(n*n) for n in range(1, K))
    if i == 4: return 1 + 2*sum((-1)**n*s**(n*n) for n in range(1, K))
    return 2*sum(exp(pi*1j*t*(n+mpf(1)/2)**2) for n in range(0, K))
def eta(t, K=400):
    q = exp(2j*pi*t); p = exp(2j*pi*t/24)
    for n in range(1, K): p *= 1 - q**n
    return p
def E2(t, K=400):
    q = exp(2j*pi*t)
    return 1 - 24*sum(sum(d for d in range(1, n+1) if n % d == 0)*q**n for n in range(1, K))
kappa = lambda t: (2*th(4, t)/th(2, t))**2 - (2*th(2, t)/th(4, t))**2
Khat = lambda t: th(3, t)**12/eta(t)**12
for tr in [mpc(0, 3), mpc(0, 6)]:
    t = 1 - 1/tr
    print("tau_r=%s  kappa=%s  Khat=%s" % (tr, mp.nstr(kappa(t), 8), mp.nstr(Khat(t), 8)))
t = mpc(0.13, 0.71); h = mpf('1e-12')
dlog = lambda f: (f(t+h) - f(t-h))/(2*h)/f(t)/(2j*pi)   # q d/dq log f
print("24 q dlog th2 - (E2+th3^4+th4^4):", mp.nstr(abs(24*dlog(lambda u: th(2, u)) - (E2(t)+th(3, t)**4+th(4, t)**4)), 5))
print("24 q dlog th3 - (E2+th2^4-th4^4):", mp.nstr(abs(24*dlog(lambda u: th(3, u)) - (E2(t)+th(2, t)**4-th(4, t)**4)), 5))
print("24 q dlog th4 - (E2-th2^4-th3^4):", mp.nstr(abs(24*dlog(lambda u: th(4, u)) - (E2(t)-th(2, t)**4-th(3, t)**4)), 5))
dk = (kappa(t+h) - kappa(t-h))/(2*h)/(2j*pi)
print("q dkappa/dq + th3^8/(th2^2 th4^2):", mp.nstr(abs(dk + th(3, t)**8/(th(2, t)**2*th(4, t)**2)), 5))
print("kappa^2 - (Khat - 64):", mp.nstr(abs(kappa(t)**2 - (Khat(t) - 64)), 5))
