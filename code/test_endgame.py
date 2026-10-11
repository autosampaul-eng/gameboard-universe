import numpy as np, crn, build, bigglass as B, sys
from scipy.spatial import cKDTree
ncell = int(sys.argv[1]); drag = int(sys.argv[2]); a = 4 * crn.D / np.sqrt(3); Lb = ncell * a; n = 8 * ncell ** 3
pos, nb = B.random_network_lean(n, Lb, 5, passes=200); deg = np.full(n, 4, np.int64); nb = nb.copy()
cand = cKDTree(pos, boxsize=Lb).query(pos, k=25)[1][:, 1:].astype(np.int32)
rng = np.random.default_rng(1); stamp = np.zeros(n, np.int64); buf = np.empty(1024, np.int64)
a0 = 0; b0 = nb[a0, 0]; build.rem(nb, deg, a0, b0); u = b0
for s in range(drag):
    x = cand[u, rng.integers(24)]
    if x == u or build.adj(nb, deg, u, x): continue
    y = nb[x, rng.integers(deg[x])]
    if y == u or y == a0: continue
    build.rem(nb, deg, x, y)
    if build.ring_le4(nb, deg, u, x, stamp, s + 1, buf): build.add(nb, deg, x, y); continue
    build.add(nb, deg, u, x); u = y
print('holes', np.nonzero(deg < 4)[0], 'separation %.1f A (box %.1f)' % (np.sqrt(B._dist2(pos, a0, u, Lb)), Lb))
steps = B.guided_endgame(nb, deg, cand, pos, Lb, 99, 10 * n)
print('guided steps', steps, 'deficient left', int((deg < 4).sum()))
if (deg == 4).all():
    sym = all(i in nb[j] for i in range(n) for j in nb[i]); dup = any(len(set(r)) < 4 for r in nb)
    bad = 0; st2 = np.zeros(n, np.int64)
    for i in range(n):
        for j in nb[i]:
            if j > i:
                build.rem(nb, deg, i, j); bad += build.ring_le4(nb, deg, i, j, st2, i * 4 + 1, buf); build.add(nb, deg, i, j)
    print('symmetric', sym, 'duplicate links', dup, 'links closing a ring <= 4:', bad)
