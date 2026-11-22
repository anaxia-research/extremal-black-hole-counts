"""Integrality of the nine extremal elliptic genera through q^6."""
from cuspcheck import genus
for m in [1, 2, 3, 4, 5, 7, 8, 11, 13]:
    phi = genus(m, 12)
    ok = all(getattr(c, 'q', 1) == 1 for k, v in phi.items() if k <= 12 for c in v.values())
    print("m=%d coefficients through q^6 integral: %s" % (m, ok), flush=True)
