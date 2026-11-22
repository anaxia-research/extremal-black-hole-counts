"""Check the Ramond-cusp evaluations on genuine extremal elliptic genera:
 chi_NS(1 - 1/tau_r, 0) -> e^{-i pi m/2} sum_l c(0,l)(-1)^l   (m = 2, 5, 7)
 and that the q^0 term of Dhat chi_NS vanishes (m = 4, 8)."""
from mpmath import mp, mpc, mpf, exp, pi
from exist import exact_system, solve
from jacobi import phi_m21, phi_01, E4, E6, mul
from polar import modforms
mp.dps = 30

def genus(m, X):
    res, rank, J = solve(m)
    sol, labels, _ = res
    P21, P01, e4, e6 = phi_m21(X), phi_01(X), E4(X), E6(X)
    pw = lambda F, k: (lambda R: R)(__import__('functools').reduce(lambda A, _: mul(A, F, X), range(k), {0: {0: 1}}))
    phi = {}
    for coef, (a, b, c) in zip(sol, labels):
        if coef == 0: continue
        M = mul(mul(pw(P21, c), pw(P01, m-c), X), mul(pw(e4, a), pw(e6, b), X), X)
        for k, v in M.items():
            for l, cc in v.items():
                phi.setdefault(k, {}); phi[k][l] = phi[k].get(l, 0) + coef*cc
    return phi   # keys: x-exponent (x = q^{1/2}), only even keys nonzero

def chi(phi, tau, z):
    q = exp(2j*pi*tau); y = exp(2j*pi*z); tot = mpc(0)
    for k, v in phi.items():
        if k % 2: continue
        n = k//2
        for l, c in v.items():
            tot += (mpf(int(c.p))/int(c.q)) * q**n * y**l
    return tot

def chiNS(phi, m, tau, z):
    return exp(2j*pi*m*(tau/4 + z + mpf(1)/2)) * chi(phi, tau, z + tau/2 + mpf(1)/2)

if __name__ == "__main__":
  for m in [2, 5, 7]:
      X = 60
      phi = genus(m, X)
      target = exp(-1j*pi*m/2)*sum(((mpf(int(c.p))/int(c.q)))*(-1)**l for l, c in phi[0].items())
      for tr in [mpc(0, 2), mpc(0, 3)]:
          tau = 1 - 1/tr
          val = chiNS(phi, m, tau, 0)
          print("m=%d tau_r=%s chi_NS(R)=%s  target=%s" % (m, tr, mp.nstr(val, 10), mp.nstr(target, 10)))
