from mpmath import mp, mpf, sqrt, cosh, exp, pi, log
mp.dps = 400
from ramond import R_explicit
from superext import beta_list
K,L=40,160
beta=beta_list(K+2)
R,e=R_explicit(K,L,beta)
def MR(k,n):
    mu=k-mpf(1)/8; mup=k-mpf(9)/8
    m=cosh(2*pi*sqrt(2*mu*n))
    if k>=2: m+=cosh(2*pi*sqrt(2*mup*n))
    return m/sqrt(n)
worst=(0,None); worst2=(0,None)
for k in range(1,K+1):
    for n in range(1,L+1):
        rel=abs(R[(k,n)]/MR(k,n)-1)
        mu=k-mpf(1)/8
        sc=rel*exp(pi*sqrt(2*mu*n)/2)   # normalised
        if rel>worst[0]: worst=(rel,(k,n))
        if sc>worst2[0]: worst2=(sc,(k,n))
print("max rel error", float(worst[0]), worst[1])
print("max rel*exp(pi sqrt(2 mu n)/2)", float(worst2[0]), worst2[1])
for (k,n) in [(1,1),(2,1),(10,1),(20,1),(40,1),(40,2),(40,5),(40,10),(40,40),(10,10)]:
    print(k,n, "%.6e"%float(R[(k,n)]/MR(k,n)-1))
w=(0,None)
for k in range(1,K+1):
    for n in range(1,L+1):
        mu=k-mpf(1)/8
        sc=abs(R[(k,n)]/MR(k,n)-1)*exp(mpf(4)/3*pi*sqrt(2*mu*n))
        if sc>w[0]: w=(sc,(k,n))
print("max rel*exp(4/3 pi sqrt(2 mu n))", float(w[0]), w[1])
import collections
