"""Random 4-regular starting network from random points (no crystal memory)."""
import numpy as np
from numba import njit
from scipy.spatial import cKDTree
import crn


@njit(cache=True)
def adj(nb, deg, x, y):
    for a in range(deg[x]):
        if nb[x, a] == y:
            return True
    return False


@njit(cache=True)
def ring_le4(nb, deg, x, y, stamp, tag, buf):
    stamp[x] = tag; buf[0] = x; n = 1; lo = 0
    for d in range(3):
        hi = n
        for q in range(lo, hi):
            u = buf[q]
            for a in range(deg[u]):
                w = nb[u, a]
                if w == y:
                    return True
                if stamp[w] != tag:
                    stamp[w] = tag; buf[n] = w; n += 1
        lo = hi
    return False


@njit(cache=True)
def add(nb, deg, x, y):
    nb[x, deg[x]] = y; deg[x] += 1; nb[y, deg[y]] = x; deg[y] += 1


@njit(cache=True)
def rem(nb, deg, x, y):
    for a in range(deg[x]):
        if nb[x, a] == y:
            nb[x, a] = nb[x, deg[x] - 1]; nb[x, deg[x] - 1] = -1; deg[x] -= 1; break
    for a in range(deg[y]):
        if nb[y, a] == x:
            nb[y, a] = nb[y, deg[y] - 1]; nb[y, deg[y] - 1] = -1; deg[y] -= 1; break


@njit(cache=True)
def greedy(pi, pj, n):
    nb = -np.ones((n, 4), np.int64); deg = np.zeros(n, np.int64)
    stamp = np.zeros(n, np.int64); buf = np.empty(1024, np.int64); tag = 0
    for k in range(pi.shape[0]):
        i = pi[k]; j = pj[k]
        if deg[i] < 4 and deg[j] < 4 and not adj(nb, deg, i, j):
            tag += 1
            if not ring_le4(nb, deg, i, j, stamp, tag, buf):
                add(nb, deg, i, j)
    return nb, deg


@njit(cache=True)
def repair(nb, deg, cand, seed, maxit):
    np.random.seed(seed)
    n = nb.shape[0]; K = cand.shape[1]
    stamp = np.zeros(n, np.int64); buf = np.empty(1024, np.int64); tag = 0
    it = 0
    for u0 in range(n):
        u = u0
        local = 0
        while deg[u] < 4 and it < maxit and local < 5000:
            it += 1; local += 1
            # 1) join a nearby deficient node
            done = False
            for t in range(K):
                v = cand[u, t]
                if v != u and deg[v] < 4 and not adj(nb, deg, u, v):
                    tag += 1
                    if not ring_le4(nb, deg, u, v, stamp, tag, buf):
                        add(nb, deg, u, v); done = True; break
            if done:
                continue
            # 2) steal: pick near x, its neighbour y; remove x-y, add u-x (y becomes deficient)
            x = cand[u, np.random.randint(min(K, 12 + local // 50))]
            if x == u or adj(nb, deg, u, x) or deg[x] == 0:
                continue
            y = nb[x, np.random.randint(deg[x])]
            if y == u:
                continue
            rem(nb, deg, x, y)
            tag += 1
            if ring_le4(nb, deg, u, x, stamp, tag, buf):
                add(nb, deg, x, y); continue
            add(nb, deg, u, x)
            tag += 1
            if deg[u] < 4 and not adj(nb, deg, u, y) and not ring_le4(nb, deg, u, y, stamp, tag, buf):
                add(nb, deg, u, y)
            else:
                u = y
    return it


def random_network(n, Lb, seed):
    rng = np.random.default_rng(seed)
    # random sequential addition with a minimum separation
    dmin = 0.62 * crn.D
    pts = np.empty((0, 3))
    while len(pts) < n:
        m = int((n - len(pts)) * 1.3) + 100
        newp = rng.uniform(0, Lb, (m, 3))
        if len(pts):
            t = cKDTree(pts, boxsize=Lb)
            d, _ = t.query(newp)
            newp = newp[d > dmin]
        t2 = cKDTree(newp, boxsize=Lb)
        keep = np.ones(len(newp), bool)
        for i, j in t2.query_pairs(dmin, output_type='ndarray'):
            if keep[i] and keep[j]:
                keep[j] = False
        pts = np.vstack([pts, newp[keep]])[:n]
    t = cKDTree(pts, boxsize=Lb)
    dd, cand = t.query(pts, k=25)
    cand = cand[:, 1:].astype(np.int64); dd = dd[:, 1:]
    pi = np.repeat(np.arange(n), 24); pj = cand.ravel(); dist = dd.ravel()
    o = np.argsort(dist + rng.uniform(0, 1e-6, dist.shape))
    nb, deg = greedy(pi[o].astype(np.int64), pj[o].astype(np.int64), n)
    print('greedy: deficient nodes', int((deg < 4).sum()), flush=True)
    it = 0
    for p in range(80):
        it += repair(nb, deg, cand, seed + p, 50 * n)
        left = int((deg < 4).sum())
        print('repair pass', p, 'iterations', it, 'deficient left', left, flush=True)
        if left == 0:
            break
    assert (deg == 4).all()
    # symmetric check
    for i in range(0, n, max(1, n // 1000)):
        for j in nb[i]:
            assert i in nb[j]
    return pts, nb


def ring_counts(nb, maxlen=4):
    """count 3- and 4-rings (each ring counted once)."""
    n = len(nb); s = [set(r) for r in nb]
    tri = sum(1 for i in range(n) for j in nb[i] if j > i for k in s[i] & s[j] if k > j)
    sq = 0
    for i in range(n):
        for j in nb[i]:
            for k in nb[j]:
                if k == i: continue
                for l in nb[k]:
                    if l == j or l == i: continue
                    if i in s[l]: sq += 1
    return tri, sq // 8
