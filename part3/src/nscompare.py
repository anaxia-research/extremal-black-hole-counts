"""Compare the Ramond bound R_k(1)=2048 e_1(k-1)-beta_k with the NS bound N_k(1) of Part II, for k<=3000."""
from superext import kappa_list, dlist, beta_list
from ramond import E_series
from math import log10
K=3000
kap=kappa_list(K+2); d=dlist(K+3); beta=beta_list(K+1)
e1=[int(v) for v in E_series(1,K).coeffs()]
worst=None
for k in range(1,K+1):
    N1=sum(m*d[k-m]*kap[m] for m in range(1,k+1))-d[k+1]
    R1=2048*e1[k-1]-beta[k]
    assert R1>N1>0,(k,R1,N1)
    lr=log10(R1)-log10(N1)
    if worst is None or lr<worst[0]: worst=(lr,k)
    if k in (1,2,4,8,9,10,100,1000,3000): print(k,"R_k(1)=%.4e N_k(1)=%.4e"%(R1,N1))
print("R_k(1)>N_k(1) for all k<=%d; min log10 ratio %.3f at k=%d"%(K,worst[0],worst[1]))
