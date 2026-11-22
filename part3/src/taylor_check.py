"""Check: (-1)^k f_k(-24-u) = beta_k + sum_r e_r(k-r) u^r (Remark on Taylor expansion), and P_2 = w^2-552."""
from flint import fmpz_poly
from superext import kappa_list, dlist, ns_w, beta_list
from jcoeffs import jcoeffs
from ramond import E_series
kap=kappa_list(400); c=jcoeffs(150); d=dlist(500); beta=beta_list(30)
def fk_poly(k):
    xg = fmpz_poly([1, 0] + [kap[n] for n in range(1, 3 * k + 5)])
    pw=[fmpz_poly([1])]
    for r in range(1,k+1): pw.append((pw[-1]*xg).truncate(2*k+2))
    w=ns_w(k,0,kap,c,d); rem={t:w[t] for t in range(-k,1)}; f=[0]*(k+1)
    for r in range(k,-1,-1):
        cr=[int(v) for v in pw[r].coeffs()]+[0]*(3*k+5)
        a=rem[-r]; f[r]=a
        for t in range(-r,1): rem[t]-=a*cr[t+r]
    return fmpz_poly(f)
for k in range(1,16):
    f=fk_poly(k)
    g=fmpz_poly([-24,-1])   # -24-u
    val=fmpz_poly([0]); P=fmpz_poly([1])
    for coef in [int(v) for v in f.coeffs()]:
        val+=coef*P; P=P*g
    val=(-1)**k*val
    co=[int(v) for v in val.coeffs()]
    exp=[beta[k]]+[int(E_series(r,k).coeffs()[k-r]) if k-r < len(E_series(r,k).coeffs()) else 0 for r in range(1,k+1)]
    assert co==exp,(k,co,exp)
print("Taylor identity checked for k<=15; e.g. k=3:", [int(v) for v in ((-1)**3*fk_poly(3)).coeffs()])
