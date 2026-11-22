from ramond import E_series, Bser, powtr
from superext import beta_list
from collections import defaultdict
K,L=40,160
B=Bser(L)
Bp={r:[int(v) for v in powtr(B,r,L).coeffs()]+[0]*(L+2) for r in range(1,min(K,L)+1)}
e={r:[int(v) for v in E_series(r,K).coeffs()]+[0]*(K+2) for r in range(1,K+1)}
def h(k,n): return sum(4096**r*Bp[r][n-r]*e[r][k-r] for r in range(1,min(k,n)+1))
def v2(x):
    v=0
    while x%2==0: x//=2; v+=1
    return v
bad=[(k,n) for k in range(1,K+1) for n in range(1,L+1) if v2(h(k,n))<v2(h(1,n))]
print("v2(h_n(k))>=v2(h_n(1)) failures:",bad[:5],len(bad))
mn=defaultdict(lambda:99)
for n in range(1,L+1): mn[n%8]=min(mn[n%8],v2(h(1,n)))
print({r:mn[r]-12 for r in range(8)})
