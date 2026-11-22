from mpmath import mp, mpf, exp, pi, sqrt, log, diff
mp.dps = 50
import cardy
from cardy import CY, Theta, D4, Pi, Ys
al = mpf(3); Y = Ys
tau0 = max(sqrt(2)*(Theta*CY + Theta*D4*exp(pi*k/2-2*pi*Y*k) + (exp(-2*pi*Y) if k==1 else 0)) for k in range(1,50))
tau0 = max(tau0, Pi*CY)
print("tau0 =", tau0)
TAU = mpf(134)
assert tau0 <= TAU
def parts(mu, r):
    n = mu*r**4; ey = exp(2*pi*mu*r**2/al**2); main = exp(2*sqrt(2)*pi*mu*r**2)/sqrt(n)
    E1 = 2*sqrt(2)*sqrt(al*r+2)*exp(pi*mu*(1+1/(al*r))**2)*ey/(pi*n)/main
    E2 = 4*sqrt(2)*sqrt(al*r+1)*exp(pi*mu*(al*r+1)**2/9)*ey/main
    E3 = 4*sqrt(al*r+1)*TAU*exp(2*pi*Y*(mu+mpf(1)/8))*ey/main
    return E1, E2, E3
def Eplus(mu, r): return sum(parts(mu, r))
for (mu, r, lab) in [(mpf(39)/8, 1, "k>=5,n>=k"), (mpf(15)/8, mpf(2)**0.25, "k>=2,n>=2k"),
                     (mpf(7)/8, (mpf(32)/7)**0.25, "k=1,n>=4"), (mpf(79)/8, 1, "k>=10,n>=k"),
                     (mpf(399)/8,1,"k>=50,n>=k"), (mpf(15)/8, 1, "k>=2,n>=k")]:
    p = parts(mu, r)
    print(lab, "E1=%.3e E2=%.3e E3=%.3e  E+=%.4e" % tuple(float(x) for x in p+(sum(p),)))
# log-derivative bounds at the corner (7/8,1), formulas as in the paper
mu0, r0 = mpf(7)/8, mpf(1)
c = 2*sqrt(2) - 2/al**2
d = [al/(2*(al*r0+2)) - 2/r0 - 2*pi*mu0*c*r0*0,  # placeholder
     ]
dr1 = al/(2*(al*r0+2)) - 2/r0 - 2*c*pi*mu0*r0*2/2
dr2 = 2/r0 + al/(2*(al*r0+1)) + pi*mu0*(2*al*(al*r0+1)/9 - 2*c*r0)
dr3 = 2/r0 + al/(2*(al*r0+1)) - 2*c*pi*mu0*r0
dm1 = pi*((1+1/(al*r0))**2 - c*r0**2)            # (the -1/(2mu) dropped)
dm2 = 1/(2*mu0) + pi*((al*r0+1)**2/9 - c*r0**2)
dm3 = 1/(2*mu0) + 2*pi*Y - pi*c*r0**2
print("corner log-derivative bounds:", [float(x) for x in (dr1, dr2, dr3, dm1, dm2, dm3)])
# grid check of monotonicity
import itertools
ok = True
for mu in [mpf(7)/8 + i*mpf(1)/2 for i in range(0, 60)]:
    prev = None
    for j in range(0, 60):
        r = 1 + mpf(j)/10
        v = Eplus(mu, r)
        if prev is not None and v > prev: ok = False
        prev = v
for r in [1 + mpf(j)/10 for j in range(0, 60)]:
    prev = None
    for i in range(0, 60):
        v = Eplus(mpf(7)/8 + i*mpf(1)/2, r)
        if prev is not None and v > prev: ok = False
        prev = v
print("E+ decreasing on grid:", ok)
for k in [5,10,20,50,100,200]:
    mu=k-mpf(1)/8; print("k=%d  E+(mu,1)=%.3e"%(k, float(Eplus(mu,1))))
