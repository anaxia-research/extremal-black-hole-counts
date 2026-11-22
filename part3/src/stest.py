from mpmath import mp, mpf, mpc, exp, pi, sqrt
mp.dps=40
from ramond import H_direct
from superext import kappa_list, dlist, ns_primaries, triangular
from jcoeffs import jcoeffs
kap=kappa_list(3000); c=jcoeffs(1500); d=dlist(3000); tri=triangular(400)
def e(z): return exp(2j*pi*z)
L=80
for k in (1,2,3):
    h=H_direct(k,L,kap,c,d)
    th=[0]*(L+1)
    for m in range(-10,11):
        if m*m<=L: th[m*m]+=(-1)**abs(m)
    Rc=[sum(th[j]*h[n-j] for j in range(n+1)) for n in range(L+1)]
    def calR(t): return sum(Rc[n]*e(n*t) for n in range(L+1))
    N,_=ns_primaries(k,60,kap,c,d,tri)
    def Pk(s):
        q=e(s); tot=e((mpf(1)/8-k)*s)-e((mpf(9)/8-k)*s)
        tot+=sum(N[j]*e((j+mpf(1)/8)*s) for j in range(1,61))
        return tot
    t0=mpc(0.1,0.9)
    lhs=calR(-1/t0)
    rhs=(-1j*t0)**0.5*sqrt(2)*(-1)**k*e(-mpf(1)/16)*Pk((t0+1)/2)
    print(k, lhs, rhs)
