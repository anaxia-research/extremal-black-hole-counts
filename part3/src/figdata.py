from math import log10
from ramond import E_series, rho_series
from superext import beta_list, kappa_list, dlist
K=40; L=40
e={r:[int(v) for v in E_series(r,K).coeffs()]+[0]*(K+2) for r in range(1,K+1)}
rho={r:[int(v) for v in rho_series(r,L).coeffs()]+[0]*(L+2) for r in range(1,L+1)}
for (k,n) in [(20,5),(20,20),(5,20)]:
    pts=[]
    for r in range(1,min(k,n)+1):
        t=4096**r*rho[r][n-r]*e[r][k-r]//2
        pts.append("(%d,%.3f)"%(r,log10(t)))
    print((k,n)," ".join(pts))
beta=beta_list(61); kap=kappa_list(70); d=dlist(70)
e1=[int(v) for v in E_series(1,61).coeffs()]
A=[];Bn=[];C=[]
for k in range(1,61):
    N1=sum(m*d[k-m]*kap[m] for m in range(1,k+1))-d[k+1]
    R1=2048*e1[k-1]-beta[k]
    if k%2==0 and beta[k]<0: A.append("(%d,%.3f)"%(k,log10(-beta[k])))
    Bn.append("(%d,%.3f)"%(k,log10(N1))); C.append("(%d,%.3f)"%(k,log10(R1)))
print("negbeta"," ".join(A)); print("N1"," ".join(Bn)); print("R1"," ".join(C))
# R_k(n) for fixed n as combos, check closed forms n=1..4
print("rho_r small:", {r:rho[r][:4] for r in range(1,5)})
