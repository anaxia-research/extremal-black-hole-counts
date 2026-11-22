import sys
from mpmath import mp, mpf, exp, pi, log, sqrt, cosh, ceil
mp.dps = 60
from ramond import R_explicit
from superext import beta_list, kappa_list
kap = kappa_list(400)
def gplus(Y): return sum(kap[N]*exp(-2*pi*N*Y) for N in range(1,401))
def gNk(Y): return sum(N*kap[N]*exp(-2*pi*N*Y) for N in range(1,401))
W = exp(2*pi)+gplus(1)
Ys = mpf('1.002')
mY = exp(2*pi*Ys)-gplus(Ys)
qY=exp(-2*pi*Ys); DY=mpf(1)
for n in range(1,60): DY*=(1+qY**(2*n+1))/(1-qY**(2*n+2))
CY = DY*(exp(2*pi*Ys)+gNk(Ys))/(mY-W)
q4=exp(-pi/2); D4=mpf(1)
for n in range(1,80): D4*=(1+q4**(2*n+1))/(1-q4**(2*n+2))
Theta = sum(exp(-pi/2*(i*(i+1)/2+mpf(1)/8)) for i in range(0,40))
xR=exp(-pi); Pi=mpf(1)
for n in range(1,200): Pi*=(1+xR**n)/(1-xR**n)
def Fk(k): return CY*exp(2*pi*k*Ys)
def Tk(k):
    TN = sqrt(2)*(Theta*(Fk(k)+exp(pi*k/2)*D4) + (1 if k==1 else 0))
    TR = Pi*Fk(k)
    return max(TN,TR)
def MR(k,n):
    mu=k-mpf(1)/8; mup=k-mpf(9)/8
    m = cosh(2*pi*sqrt(2*mu*n))
    if k>=2: m += cosh(2*pi*sqrt(2*mup*n))
    return m/sqrt(n)
def Ebound(k,n,alpha=3):
    mu=k-mpf(1)/8; mup=k-mpf(9)/8
    r=(n/mu)**(mpf(1)/4); Q=int(ceil(alpha*r)); Q=max(Q,2)
    y=mpf(1)/Q**2; ey=exp(2*pi*n*y)
    ells=[mu]+([mup] if k>=2 else [])
    R1=sum(sqrt(2)*sqrt(Q+1)*exp(pi*l*(1+mpf(1)/Q)**2) for l in ells)*ey/(pi*n)
    R2=4*sqrt(2)*(sqrt(Q)-1)*exp(pi*mu*Q**2/9)*ey if Q>=3 else 0
    R3=4*sqrt(Q)*Tk(k)*ey
    return (R1+R2+R3)/(2*MR(k,n)), Q, (R1/(2*MR(k,n)),R2/(2*MR(k,n)),R3/(2*MR(k,n)))
if __name__=="__main__":
    print("C_Y=",CY," Theta=",Theta," Pi=",Pi)
    K,L=30,240
    beta=beta_list(K+2)
    R,_=R_explicit(K,L,beta)
    for k in [1,2,3,5,10,20,30]:
        for lam in [1,2,4,8]:
            n=lam*k if lam*k>=1 else 1
            if n>L: continue
            true=abs(R[(k,n)]/MR(k,n)-1)
            E,Q,parts=Ebound(k,n)
            print(k,n,"true %.3e  E %.3e Q=%d"%(float(true),float(E),Q), ["%.1e"%float(p) for p in parts])
