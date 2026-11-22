from cuspcheck import genus
for m in [1, 2, 3, 4, 5, 7, 8, 11, 13]:
    phi = genus(m, 8)
    for n in [1, 2]:
        row = [phi.get(2*n, {}).get(l, 0) for l in range(0, m+2*n+1)]
        while row and row[-1] == 0: row.pop()
        print("m=%d n=%d" % (m, n), [str(c) for c in row], "integral:", all(getattr(c, 'q', 1) == 1 for c in row))
