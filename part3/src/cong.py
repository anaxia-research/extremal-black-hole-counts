from ramond import R_explicit, H_direct, E_series, rho_series
from superext import beta_list, kappa_list, dlist
from jcoeffs import jcoeffs
K,L=40,120
beta=beta_list(K+2)
R,e=R_explicit(K,L,beta)
def v2(x):
    if x==0: return 99
    v=0
    while x%2==0: x//=2; v+=1
    return v
# R mod 4096 vs (-1)^m beta at squares, 0 elsewhere, when n != 1 mod 8
import math
bad=[]
for k in range(1,K+1):
    for n in range(1,L+1):
        m=math.isqrt(n); sq=(m*m==n)
        base=((-1)**m*beta[k]) if sq else 0
        d=R[(k,n)]-base
        if (n-1)%8!=0 and d%4096: bad.append((k,n))
        if d%2048: bad.append(('2048',k,n))
print("R congruence mod 4096 for n!=1 mod 8 failures:",bad[:10])
# minimal 2-adic valuation of R_k(n)-base by residue of n mod 8
from collections import defaultdict
mv=defaultdict(lambda:99)
for k in range(1,K+1):
    for n in range(1,L+1):
        m=math.isqrt(n); base=((-1)**m*beta[k]) if m*m==n else 0
        mv[(n%8)]=min(mv[n%8],v2(R[(k,n)]-base))
print("min v2(R-base) by n mod 8:",dict(sorted(mv.items())))
kap=kappa_list(400); c=jcoeffs(150); d=dlist(500)
mh=defaultdict(lambda:99)
for k in range(1,16):
    h=H_direct(k,60,kap,c,d)
    for n in range(1,61): mh[n%8]=min(mh[n%8],v2(h[n]))
print("min v2(h_n) by n mod 8:",dict(sorted(mh.items())))
# parity of e_1
e1=[int(v) for v in E_series(1,200).coeffs()]
print("e_1 mod 2 first 64:", ''.join(str(x%2) for x in e1[:64]))
r1=[int(v) for v in rho_series(1,200).coeffs()]
print("rho_1 odd positions <200:", [i for i in range(200) if r1[i]%2][:30])
