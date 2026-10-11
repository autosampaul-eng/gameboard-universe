"""Memory-lean, resumable versions of build.random_network and crn.global_relax (same recipe, bit-identical results).
- random_network_lean: chunked neighbour query, int32 candidate table, greedy pass reading pairs through the sort order
  (no full-size pi/pj copies). Random-number calls are the same, in the same order, as build.random_network.
- FIRE relaxation (crn.fire with mov = cen = all atoms) run in chunks with its full state carried between calls,
  so it can be paused and resumed across jobs and gives the same trajectory as one uninterrupted run."""
import numpy as np, os, time, json
from numba import njit
from scipy.spatial import cKDTree
import crn, build

@njit(cache=True)
def greedy_lean(o, cflat, K, n):
    nb = -np.ones((n, 4), np.int64); deg = np.zeros(n, np.int64)
    stamp = np.zeros(n, np.int64); buf = np.empty(1024, np.int64); tag = 0
    for k in range(o.shape[0]):
        q = o[k]; i = q // K; j = cflat[q]
        if deg[i] < 4 and deg[j] < 4 and not build.adj(nb, deg, i, j):
            tag += 1
            if not build.ring_le4(nb, deg, i, j, stamp, tag, buf):
                build.add(nb, deg, i, j)
    return nb, deg

def random_network_lean(n, Lb, seed, passes=80, chunk=1_000_000):
    rng = np.random.default_rng(seed)
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
    K = 24; cand = np.empty((n, K), np.int32); dist = np.empty(n * K)
    for s in range(0, n, chunk):
        e = min(n, s + chunk); dd, cc = t.query(pts[s:e], k=K + 1)
        cand[s:e] = cc[:, 1:]; dist[s * K:e * K] = dd[:, 1:].ravel()
    del t; dist += rng.uniform(0, 1e-6, dist.shape)
    o = np.argsort(dist); del dist
    nb, deg = greedy_lean(o, cand.ravel(), K, n); del o
    print('greedy: deficient nodes', int((deg < 4).sum()), flush=True)
    it = 0
    cand64 = cand  # repair reads cand[u, t]; values identical (int32)
    for p in range(passes):
        it += build.repair(nb, deg, cand64, seed + p, 50 * n)
        left = int((deg < 4).sum())
        print('repair pass', p, 'iterations', it, 'deficient left', left, flush=True)
        if left == 0:
            break
    left = int((deg < 4).sum())
    if left:
        steps = guided_endgame(nb, deg, cand, pts, Lb, seed + 10_000, 10 * n)
        print('guided endgame: deficient', left, '-> steps', steps, 'deficient left', int((deg < 4).sum()), flush=True)
    assert (deg == 4).all()
    for i in range(0, n, max(1, n // 1000)):
        for j in nb[i]:
            assert i in nb[j]
    return pts, nb

@njit(cache=True)
def fire_chunk(pos, nb, Lb, g, v, st, ftol, maxit, nits, dmax=1e300):
    """Continue crn.fire (mov = cen = all atoms, Ecut off) for up to nits iterations.
    st = [dt, alpha, npos, it, E, fmax, converged, clamped atom-steps]. dmax = largest move per atom per step."""
    n = pos.shape[0]; dtmax = 0.3
    dt = st[0]; alpha = st[1]; npos = int(st[2]); it0 = int(st[3]); E = st[4]; fmax = st[5]
    for it in range(it0, min(maxit, it0 + nits)):
        for x in range(n):
            for c in range(3):
                g[x, c] = 0.0
            for a in range(4):
                y = nb[x, a]
                for c in range(3):
                    g[y, c] = 0.0
        E = 0.0
        for x in range(n):
            E += crn.atom_eg(x, pos, nb, Lb, g)
        P = 0.0; vn = 0.0; fn = 0.0; fmax = 0.0
        for x in range(n):
            for c in range(3):
                f = -g[x, c]
                P += f * v[x, c]; vn += v[x, c] ** 2; fn += f * f
                if abs(f) > fmax:
                    fmax = abs(f)
        if fmax < ftol:
            st[0] = dt; st[1] = alpha; st[2] = npos; st[3] = it; st[4] = E; st[5] = fmax; st[6] = 1.0
            return
        vn = np.sqrt(vn); fn = np.sqrt(fn) + 1e-30
        if P > 0:
            for x in range(n):
                for c in range(3):
                    v[x, c] = (1 - alpha) * v[x, c] + alpha * vn * (-g[x, c]) / fn
            npos += 1
            if npos > 5:
                dt = min(dt * 1.1, dtmax); alpha *= 0.99
        else:
            npos = 0; dt *= 0.5; alpha = 0.1
            for x in range(n):
                for c in range(3):
                    v[x, c] = 0.0
        for x in range(n):
            s2 = 0.0
            for c in range(3):
                v[x, c] += -g[x, c] * dt
                s2 += (v[x, c] * dt) ** 2
            if s2 > dmax * dmax:          # safeguard: no atom moves more than dmax in one step (counted)
                sc = dmax / np.sqrt(s2); st[7] += 1.0
                for c in range(3):
                    pos[x, c] += v[x, c] * dt * sc
            else:
                for c in range(3):
                    pos[x, c] += v[x, c] * dt
    st[0] = dt; st[1] = alpha; st[2] = npos; st[3] = min(maxit, it0 + nits); st[4] = E; st[5] = fmax
    if st[3] >= maxit: st[6] = -1.0          # hit the iteration cap (as crn.fire would)

def new_state():
    return np.array([0.05, 0.1, 0.0, 0.0, 0.0, np.inf, 0.0, 0.0])

def relax(pos, nb, Lb, v, st, ftol=1e-3, maxit=5000, budget_s=1e30, chunk=20, log=print, dmax=1e300):
    g = np.zeros_like(pos); t0 = time.time()
    while st[6] == 0 and time.time() - t0 < budget_s:
        fire_chunk(pos, nb, Lb, g, v, st, ftol, maxit, chunk, dmax)
        if int(st[3]) % 200 == 0 or st[6] != 0: log(f"FIRE it {int(st[3])} E/atom {st[4]/len(pos):.5f} fmax {st[5]:.2e} dt {st[0]:.3f} clamped {int(st[7])} ({time.time()-t0:.0f}s)")
    return st

@njit(cache=True)
def _mi(d, Lb):
    return d - Lb * np.round(d / Lb)

@njit(cache=True)
def _dist2(pos, i, j, Lb):
    s = 0.0
    for c in range(3):
        d = _mi(pos[j, c] - pos[i, c], Lb); s += d * d
    return s

@njit(cache=True)
def guided_endgame(nb, deg, cand, pos, Lb, seed, maxsteps):
    """Pair up the remaining deficient nodes: walk the deficiency of u toward its nearest deficient partner t with the
    repair's own steal move (remove x-y, add u-x, deficiency passes to y), choosing x among u's candidates closest to t.
    Same constraints as build.repair: no duplicate links, no rings of length <= 4."""
    np.random.seed(seed)
    n = nb.shape[0]; K = cand.shape[1]
    stamp = np.zeros(n, np.int64); buf = np.empty(1024, np.int64); tag = 0
    steps = 0
    while steps < maxsteps:
        holes = np.nonzero(deg < 4)[0]
        if holes.shape[0] == 0:
            return steps
        u = holes[0]; t = holes[1] if holes.shape[0] > 1 else holes[0]
        best = 1e300
        for h in holes[1:]:
            d2 = _dist2(pos, u, h, Lb)
            if d2 < best: best = d2; t = h
        while deg[u] < 4 and steps < maxsteps:
            steps += 1
            done = False
            for q in range(K):                       # 1) join a nearby deficient node (as in repair)
                v = cand[u, q]
                if v != u and deg[v] < 4 and not build.adj(nb, deg, u, v):
                    tag += 1
                    if not build.ring_le4(nb, deg, u, v, stamp, tag, buf):
                        build.add(nb, deg, u, v); done = True; break
            if done:
                break
            # 2) guided steal: x among the 3 candidates of u closest to t (random pick), y a random neighbour of x
            #    half the moves are the repair's own random steal, so the walk can step around blocked joins
            if np.random.random() < 0.5:
                x = cand[u, np.random.randint(12)]
            else:
                dd = np.empty(K)
                for q in range(K): dd[q] = _dist2(pos, cand[u, q], t, Lb)
                order = np.argsort(dd)
                x = cand[u, order[np.random.randint(3)]]
            if x == u or build.adj(nb, deg, u, x) or deg[x] == 0:
                continue
            y = nb[x, np.random.randint(deg[x])]
            if y == u:
                continue
            build.rem(nb, deg, x, y)
            tag += 1
            if build.ring_le4(nb, deg, u, x, stamp, tag, buf):
                build.add(nb, deg, x, y); continue
            build.add(nb, deg, u, x)
            u = y
    return -steps
