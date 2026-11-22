from mpmath import mp, mpf, exp, pi, log, sqrt, nsum, inf
mp.dps = 40
from superext import kappa_list
kap = kappa_list(400)
def gplus(Y): return sum(kap[N]*exp(-2*pi*N*Y) for N in range(1,401))
def gNk(Y): return sum(N*kap[N]*exp(-2*pi*N*Y) for N in range(1,401))
gi = exp(2*pi) + gplus(1)
xR = exp(-pi)
B = mpf(1)
for n in range(1,200): B *= (1+xR**n)**24
WR = 24 + 4096*xR*B
Pi = mpf(1)
for n in range(1,200): Pi *= (1+xR**n)/(1-xR**n)
Theta = sum(exp(-pi/2*(i*(i+1)/2+mpf(1)/8)) for i in range(0,40))
print("g(i)=",gi," W_R=",WR," Pi=",Pi," Theta_1/4=",Theta)
W = max(gi, WR)
for Y in [mpf('1.0005'),mpf('1.001'),mpf('1.002'),mpf('1.005')]:
    mY = exp(2*pi*Y)-gplus(Y)
    D = mpf(1)
    # D(q)=prod (1+q^{2n+1})/(1-q^{2n+2})
    q=exp(-2*pi*Y)
    for n in range(1,50): D*= (1+q**(2*n+1))/(1-q**(2*n+2))
    C = D*(exp(2*pi*Y)+gNk(Y))/(mY-W)
    print(Y, "m_Y=",mY, "m_Y-W=",mY-W, "C=",C)
