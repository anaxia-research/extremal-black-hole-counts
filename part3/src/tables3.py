from mpmath import mp, mpf, sqrt, cosh, exp, pi, log10, log
mp.dps = 400
from ramond import R_explicit, E_series, rho_series
from superext import beta_list
import envelope
mp.dps = 400
K,L=40,160
beta=beta_list(K+2)
R,e=R_explicit(K,L,beta)
print("=== R_k(n), k<=12, n<=5")
R12,_=R_explicit(12,6,beta)
for k in range(1,13): print(k, beta[k], [R12[(k,n)] for n in range(1,6)])
print("=== e_r(N)")
for r in (1,2,3): print(r, [int(v) for v in E_series(r,10).coeffs()][:11])
print("=== rho_r(m)")
for r in (1,2): print(r, [int(v) for v in rho_series(r,12).coeffs()][:13])
def MR(k,n):
    mu=k-mpf(1)/8; mup=k-mpf(9)/8
    m=cosh(2*pi*sqrt(2*mu*n))
    if k>=2: m+=cosh(2*pi*sqrt(2*mup*n))
    return m/sqrt(n)
print("=== accuracy table")
for (k,n) in [(1,1),(1,4),(1,16),(2,1),(2,4),(2,16),(5,5),(5,20),(10,10),(10,40),(20,20),(40,40),(40,160)]:
    rel=R[(k,n)]/MR(k,n)-1
    mu=k-mpf(1)/8; r=(n/mu)**(mpf(1)/4)
    Ep = envelope.Eplus(mu, r) if r>=1 else None
    print(k,n,"%.4e"%float(rel), ("%.3e"%float(Ep)) if Ep is not None else "-")
print("=== plot data: sqrt(2 mu n), log10|rel|")
for k in (1,2,5,10,20,40):
    pts=[]
    for n in list(range(1,11))+list(range(12,161,8)):
        mu=k-mpf(1)/8
        rel=abs(R[(k,n)]/MR(k,n)-1)
        pts.append("(%.3f,%.3f)"%(float(sqrt(2*mu*n)), float(log10(rel))))
    print(k, " ".join(pts))
