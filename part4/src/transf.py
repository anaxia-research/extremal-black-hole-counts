"""Numerical check of the transformation laws used on the major arc (s = e^{pi i tau})."""
from mpmath import mp, mpc, exp, pi, mpf, tan, sqrt
mp.dps = 40
K = 400
def s_of(t): return exp(pi*1j*t)
def G(t):
    s = s_of(t); p = mpf(1)
    for n in range(1, K):
        p *= 1/((1 - s**(2*n-1))*(1 - s**(4*n)))**2
    return p
def th3(t):
    s = s_of(t); return 1 + 2*sum(s**(n*n) for n in range(1, 60))
def th4(t):
    s = s_of(t); return 1 + 2*sum((-1)**n*s**(n*n) for n in range(1, 60))
def th2(t):
    return 2*sum(exp(pi*1j*t*(n+mpf(1)/2)**2) for n in range(0, 60))
def T(t): return th4(t)**4 - th2(t)**4
def U(t):
    s = s_of(t); return sum(s**(2*n+1)/(1+s**(2*n+1))**2 for n in range(0, K))
def P3(t):
    s = s_of(t); p = mpf(1)
    for n in range(1, K): p *= (1 - s**(2*n))**3
    return p
for t in [mpc(0.03, 0.05), mpc(-0.02, 0.04), mpc(0.1, 0.3)]:
    tp = -1/t; w = -1j*t
    g1 = G(t); g2 = exp(pi*1j*t/4)*w*exp(pi*1j/(4*t))*G(tp)
    t1 = T(t); t2 = -(w**-2)*T(tp)
    u1 = U(t); u2 = -1/(4*pi*1j*t) + t**-2*U(tp)
    gp = th3(t)/P3(t)
    print("G:", mp.nstr(abs(g1/g2-1), 5), " G=th3/P3:", mp.nstr(abs(gp/g1-1), 5), " T:", mp.nstr(abs(t1/t2-1), 5), " U:", mp.nstr(abs(u1/u2-1), 5))
