"""Constants of the a priori bounds (Lemma on a priori bounds) and of the envelope E^+ (Theorem E)."""
from mpmath import mp, mpf, exp, pi, sqrt, log
mp.dps = 50
from superext import kappa_list
kap = kappa_list(600)
def gplus(Y): return sum(kap[N]*exp(-2*pi*N*Y) for N in range(1,601))
def gN(Y): return sum(N*kap[N]*exp(-2*pi*N*Y) for N in range(1,601))
# tail of g^+ beyond 600 at Y>=1/4 is below e^{-2 pi 600 (1/4)} * g(i/8)-ish: negligible; we bound it by kappa(N)<=e^{4 pi sqrt N}
tail = sum(exp(4*pi*sqrt(N)-2*pi*N*mpf(1)/4) for N in range(601,3000))
gi = exp(2*pi)+gplus(1)
xR = exp(-pi); Bv = mpf(1)
for n in range(1,300): Bv *= (1+xR**n)**24
WR = 24+4096*xR*Bv
W = mpf('536.015')
print("g(i) = %s  W_R = %s  W = %s  tail<%s" % (mp.nstr(gi,12), mp.nstr(WR,12), W, mp.nstr(tail,3)))
assert gi < W and WR < W
Y = mpf('1.002')
mY = exp(2*pi*Y)-gplus(Y)
q=exp(-2*pi*Y); DY=mpf(1)
for n in range(1,80): DY*=(1+q**(2*n+1))/(1-q**(2*n+2))
C = DY*(exp(2*pi*Y)+gN(Y))/(mY-W)
print("m_Y = %s, D(e^{-2 pi Y}) = %s, C_* = %s" % (mp.nstr(mY,10), mp.nstr(DY,12), mp.nstr(C,10)))
q4=exp(-pi/2); D4=mpf(1)
for n in range(1,200): D4*=(1+q4**(2*n+1))/(1-q4**(2*n+2))
Theta = sum(exp(-pi/2*(mpf(i*(i+1))/2+mpf(1)/8)) for i in range(60))
Pi = mpf(1)
for n in range(1,300): Pi*=(1+xR**n)/(1-xR**n)
print("Theta = %s  D(e^{-pi/2}) = %s  Pi = %s" % (mp.nstr(Theta,12), mp.nstr(D4,10), mp.nstr(Pi,10)))
T = max(max(sqrt(2)*(Theta*(C+D4*exp(pi*k/2-2*pi*Y*k))+(sqrt(2)*exp(-2*pi*Y) if k==1 else 0)) for k in range(1,200)), Pi*C)
print("tau = max T_k e^{-2 pi Y k} = %s (<= 134)" % mp.nstr(T,10)); assert T <= 134
# envelope
al=mpf(3); c=2*sqrt(2)-2/al**2
def parts(mu,r):
    E1=2*sqrt(2)/pi*sqrt(al*r+2)/(sqrt(mu)*r**2)*exp(pi*mu*((1+1/(al*r))**2-c*r**2))
    E2=4*sqrt(2)*sqrt(al*r+1)*sqrt(mu)*r**2*exp(pi*mu*((al*r+1)**2/9-c*r**2))
    E3=4*134*sqrt(al*r+1)*sqrt(mu)*r**2*exp(2*pi*Y*(mu+mpf(1)/8)-c*pi*mu*r**2)
    return E1,E2,E3
def Ep(mu,r): return sum(parts(mu,r))
print("c =", mp.nstr(c,10))
for lab,mu,r in [("k>=5",mpf(39)/8,1),("k>=10",mpf(79)/8,1),("k>=50",mpf(399)/8,1),("k>=2,n>=2k",mpf(15)/8,mpf(2)**0.25),("k=1,n>=4",mpf(7)/8,(mpf(32)/7)**0.25)]:
    p=parts(mu,r); print(lab, " ".join(mp.nstr(x,5) for x in p), " E+ =", mp.nstr(sum(p),8))
mu0=mpf(7)/8; r0=mpf(1)
vals=[al/(2*(al*r0+2))-2/r0-2*c*pi*mu0*r0, 2/r0+al/(2*(al*r0+1))+pi*mu0*(2*al*(al*r0+1)/9-2*c*r0), 2/r0+al/(2*(al*r0+1))-2*c*pi*mu0*r0,
      pi*((1+1/(al*r0))**2-c*r0**2), 1/(2*mu0)+pi*((al*r0+1)**2/9-c*r0**2), 1/(2*mu0)+2*pi*Y-pi*c*r0**2]
print("corner log-derivative bounds:", [mp.nstr(v,5) for v in vals]); assert all(v<0 for v in vals)
worst=max(Ep(mu,1)/(sqrt(mu)*exp(-mpf('0.6')*pi*mu)) for mu in [mpf(7)/8+mpf(j)/8 for j in range(0,4000)])
print("max E+(mu,1)/(sqrt(mu) e^{-0.6 pi mu}) over 7/8<=mu<=500: %s (<=2400)" % mp.nstr(worst,8)); assert worst<=2400
