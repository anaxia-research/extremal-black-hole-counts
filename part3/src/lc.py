from ramond import R_explicit, E_series, Bser, powtr
from superext import beta_list
from collections import defaultdict
K,L=40,160
beta=beta_list(K+2)
R,e=R_explicit(K,L,beta)
bad=[(k,n) for k in range(1,K+1) for n in range(2,L) if not R[(k,n)]**2>R[(k,n-1)]*R[(k,n+1)]]
print("log-concavity in n failures:",bad[:10], len(bad))
badk=[(k,n) for k in range(2,K) for n in range(1,L+1) if not R[(k,n)]**2>R[(k-1,n)]*R[(k+1,n)]]
print("log-concavity in k failures:",badk[:10], len(badk))
# h_n via explicit formula
B=Bser(L)
Bp={r:[int(v) for v in powtr(B,r,L).coeffs()]+[0]*(L+2) for r in range(1,min(K,L)+1)}
def h(k,n): return sum(4096**r*Bp[r][n-r]*e[r][k-r] for r in range(1,min(k,n)+1))
def v2(x):
    v=0
    while x%2==0: x//=2; v+=1
    return v
mv=defaultdict(lambda:99); arg={}
for k in range(1,K+1):
    for n in range(1,L+1):
        x=h(k,n); assert x>0
        if v2(x)<mv[n%8]: mv[n%8]=v2(x); arg[n%8]=(k,n)
print("min v2(h_n) by n mod 8 (k<=40,n<=160):",dict(sorted(mv.items())), arg)
badlc=[(k,n) for k in range(1,K+1) for n in range(2,L) if not h(k,n)**2>h(k,n-1)*h(k,n+1)]
print("h log-concave failures:", badlc[:5], len(badlc))
