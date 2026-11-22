"""Effective circle-method bounds for the even cases (m = 2 mod 4: L(m) < -12m;
m = 0 mod 4: R1 > 0, R2 > 0), valid for N = m/2 >= N0.  All constants are
evaluated with 50 digits and rounded in the safe direction by margins."""
from mpmath import mp, mpf, sqrt, exp, pi, log, cos, sin, erf, factorial, binomial, quad
mp.dps = 50
from flint import fmpz_poly
from sseries import series, inv_series

# ---------- constants of the q'-expansions at Im tau' >= 1/2 ----------
r0 = exp(-pi/2)
S = 400
G, T, D = series(S)
GT = G.mul_low(T, S+1)
def sigma(n): return sum(d for d in range(1, n+1) if n % d == 0)
def absum(poly, start):   # sum_{n>=start} |c_n| r0^(n-start)
    return sum(abs(int(poly[n]))*r0**(n-start) for n in range(start, S+1))
C_G = absum(G, 1); G_r0 = 1 + C_G*r0
C_GT = absum(GT, 1); GT_max = abs(int(GT[0])) + C_GT*r0
C_U = sum(sigma(j)*r0**(j-1) for j in range(1, S+1))
# tails beyond S are < 1e-40 (coefficients grow like e^{pi sqrt n}, r0^n = e^{-pi n/2})
print("C_G=%s G(r0)=%s C_GT=%s GT_max=%s C_U=%s" % tuple(mp.nstr(v, 8) for v in (C_G, G_r0, C_GT, GT_max, C_U)))

# ---------- amplitudes alpha(w), w = -i tau: Taylor data and Cauchy bounds ----------
# alpha_L  = -tanh(pi w/2)/w * e^{-pi w/4},     j0 = 0, a0 = -pi/2
# alpha_R2 = (1/(4 pi)) tanh(pi w/2) e^{-pi w/4}, j0 = 1, a1 = 1/8
# alpha_R1 = (1/2) w tanh(pi w/2) sech^2(pi w/2) e^{-pi w/4}, j0 = 2, a2 = pi/4
# On |w| = 1/2: |tanh| <= 1, |sech| <= sqrt 2, |e^{-pi w/4}| <= e^{pi/8}.
Rp = mpf(1)/2
M_L = (1/Rp)*exp(pi/8); M_R2 = exp(pi/8)/(4*pi); M_R1 = Rp*mpf(1)/2*2*exp(pi/8)
amps = {
  'L':  dict(j0=0, a=-pi/2,  M=M_L,  H_C=C_GT, H_max=GT_max, extra=False),
  'R2': dict(j0=1, a=mpf(1)/8, M=M_R2, H_C=C_G,  H_max=G_r0,  extra=True),
  'R1': dict(j0=2, a=pi/4,   M=M_R1, H_C=C_G,  H_max=G_r0,  extra=False),
}

def In_lower(n, z):
    th0 = z**(-mpf(1)/4)
    return exp(z)/sqrt(2*pi*z)*cos(n*th0)*erf(z**(mpf(1)/4)/sqrt(2)) - exp(z*cos(th0))

def bounds(name, N):
    A = amps[name]; N = mpf(N)
    eps = 1/(2*sqrt(N)); delta = sqrt(2*eps - eps**2); z = pi*sqrt(N); r = exp(-pi*eps)
    j0, J = A['j0'], A['j0']+1
    main = abs(A['a'])*pi*eps**(j0+1)*In_lower(j0+1, z)
    # ray error of S_{j0}
    ray = exp(pi/8)*exp(pi*N*eps)*(eps*(eps+delta)**j0 + sum(binomial(j0, i)*delta**(j0-i)*factorial(i)/(pi*N)**(i+1) for i in range(j0+1)))
    ray *= abs(A['a'])
    # Taylor remainder: |alpha - P_J| <= M (|w|/R')^J/(1-|w|/R'), |w| <= sqrt(2 eps)
    om = sqrt(2*eps)
    CJ = A['M']/Rp**J/(1 - om/Rp)
    KJ = 2**(mpf(J)/2)*eps**(J+1)*sqrt(8*eps)*exp(pi/(4*eps)) + 2*delta*(2*eps)**(mpf(J)/2)*exp(pi/(8*eps))
    ET = mpf(1)/2*CJ*A['H_max']*exp(pi*N*eps)*KJ
    # P_J (H-1) term
    PJmax = abs(A['a'])*om**j0
    EH = delta*PJmax*A['H_C']*exp(-3*pi/8)*exp(pi*N*eps)
    EX = 0
    if A['extra']:
        EX = delta*2*C_U*G_r0*exp(-3*pi/8)*exp(pi*N*eps)
    # minor arcs
    a_ = 1/r - r; b_ = 1/r + r; th = pi*delta
    Dmin = 1/a_ - a_*cos(th)/(a_**2*cos(th)**2 + b_**2*sin(th)**2)
    Gr = eps*exp(pi/(4*eps))*(1 + C_G*exp(-pi/eps))
    if name == 'L':
        th3r = eps**(-mpf(1)/2)*(1 + 3*exp(-pi/eps))
        amp = 2/(1-r)*2*th3r**4
    elif name == 'R2':
        amp = 2/(1-r)/(8*eps**2)
    else:
        amp = 4/(1-r)**3
    Emin = amp*Gr*exp(-2*Dmin)*exp(pi*N*eps)
    extra_target = 24*N if name == 'L' else 0
    err = ray + ET + EH + EX + Emin + extra_target
    return main, err, dict(ray=ray/main, ET=ET/main, EH=EH/main, EX=EX/main, Emin=Emin/main, target=extra_target/main, Dmin_eps=Dmin*eps)

if __name__ == "__main__":
    N0 = 20000
    for name in ['L', 'R2', 'R1']:
        main, err, parts = bounds(name, N0)
        print(name, "N0=%d  error/main = %s" % (N0, mp.nstr(err/main, 6)), {k: mp.nstr(v, 3) for k, v in parts.items()})
        worst = 0
        for e in range(43, 241):     # N from 2e4 to 1e24 on a log grid
            N = int(mpf(10)**(mpf(e)/10)) if e > 43 else N0
            m_, er, _ = bounds(name, max(N, N0))
            worst = max(worst, er/m_)
        print("   max error/main on log grid N in [2e4, 1e24]:", mp.nstr(worst, 6))
