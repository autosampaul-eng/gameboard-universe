"""Experiment 2: light cone + light clock on 3D networks (glass, diamond, cubic)."""
import sys, json, time
import numpy as np
from numba import njit
from scipy.spatial import cKDTree
import crn


def tile(pos, nb, Lb, T):
    n = len(pos)
    r = pos[nb] - pos[:, None, :]
    r -= Lb * np.round(r / Lb)
    unw = pos[:, None, :] + r
    shift = np.round((unw - pos[nb]) / Lb).astype(np.int64)  # image of neighbour
    # FIX (Oct 7, 2026): order cells to match cid = i + T*j + T*T*k (was i-major, which misplaced cells with i != k)
    cells = np.array([[q % T, (q // T) % T, q // (T * T)] for q in range(T ** 3)])
    cid = lambda c: c[..., 0] + T * c[..., 1] + T * T * c[..., 2]
    P = np.concatenate([pos + Lb * c for c in cells])
    NB = np.empty((n * len(cells), nb.shape[1]), np.int64)
    for q, c in enumerate(cells):
        tc = (c[None, None, :] + shift) % T
        NB[q * n:(q + 1) * n] = nb + n * cid(tc)
    return P, NB, Lb * T


def cubic(m):
    g = np.array([[i, j, k] for i in range(m) for j in range(m) for k in range(m)], float)
    idx = lambda i, j, k: (i % m) * m * m + (j % m) * m + (k % m)
    nb = np.array([[idx(i + 1, j, k), idx(i - 1, j, k), idx(i, j + 1, k), idx(i, j - 1, k), idx(i, j, k + 1), idx(i, j, k - 1)]
                   for i in range(m) for j in range(m) for k in range(m)], np.int64)
    return g, nb, float(m)


@njit(cache=True)
def bfs(nb, src, maxd, dist, buf):
    """dist must be filled with -1 on entry; returns number of visited nodes in buf (for reset)."""
    dist[src] = 0; buf[0] = src; m = 1; lo = 0
    for d in range(maxd):
        hi = m
        for q in range(lo, hi):
            u = buf[q]
            for c in range(nb.shape[1]):
                w = nb[u, c]
                if dist[w] < 0:
                    dist[w] = d + 1; buf[m] = w; m += 1
        lo = hi
        if lo == m:
            break
    return m


class Net:
    def __init__(self, pos, nb, Lb):
        self.pos, self.nb, self.Lb = pos, nb, Lb
        self.tree = cKDTree(pos % Lb, boxsize=Lb)
        self.dist = -np.ones(len(pos), np.int64); self.buf = np.empty(len(pos), np.int64)

    def node(self, x):
        return self.tree.query(np.asarray(x) % self.Lb)[1]

    def run_bfs(self, src, maxd):
        m = bfs(self.nb, src, maxd, self.dist, self.buf)
        return m

    def reset(self, m):
        self.dist[self.buf[:m]] = -1


def fib_dirs(k):
    i = np.arange(k) + 0.5
    phi = np.arccos(1 - 2 * i / k); th = np.pi * (1 + 5 ** 0.5) * i
    return np.c_[np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)]


def light_cone(net, centers, radii, ndir=200, cone_deg=12):
    dirs = fib_dirs(ndir); cosmax = np.cos(np.radians(cone_deg))
    out = {r: np.zeros(ndir) for r in radii}; cnt = {r: np.zeros(ndir) for r in radii}
    rmax = max(radii) + 6
    for c in centers:
        s = net.node(c)
        m = net.run_bfs(s, 10 ** 6 if False else 400)
        vis = net.buf[:m]; d = net.dist[vis]
        rel = net.pos[vis] - net.pos[s]; rel -= net.Lb * np.round(rel / net.Lb)
        rr = np.linalg.norm(rel, axis=1)
        for r in radii:
            sel = (rr > r - 3) & (rr < r + 3) & (d > 0)
            u = rel[sel] / rr[sel, None]; sp = rr[sel] / d[sel]
            dots = u @ dirs.T
            for k in range(ndir):
                w = dots[:, k] > cosmax
                if w.any():
                    out[r][k] += sp[w].mean(); cnt[r][k] += 1
        net.reset(m)
    res = {}
    for r in radii:
        cdir = out[r] / np.maximum(cnt[r], 1)
        res[r] = dict(c_mean=float(cdir.mean()), aniso=float((cdir.max() - cdir.min()) / cdir.mean()), c_dir=cdir.tolist())
    return res, dirs


def c_along(net, u, centers, r=70.0):
    """speed limit along direction u: straight-line distance / steps for nodes near the ray at distance r."""
    sp = []
    for c in centers:
        s = net.node(c)
        m = net.run_bfs(s, 400)
        vis = net.buf[:m]; d = net.dist[vis]
        rel = net.pos[vis] - net.pos[s]; rel -= net.Lb * np.round(rel / net.Lb)
        rr = np.linalg.norm(rel, axis=1)
        sel = (rr > r - 4) & (rr < r + 4) & (d > 0)
        cos = (rel[sel] @ u) / rr[sel]
        w = cos > np.cos(np.radians(6))
        if w.any():
            sp.append((rr[sel][w] / d[sel][w]).mean())
        net.reset(m)
    return float(np.mean(sp))


def tick(net, pA, pB, vel, tmax=600):
    a0 = net.node(pA)
    m = net.run_bfs(a0, tmax)
    t = 1; hit = -1
    while t <= tmax:
        nbB = net.node(pB + vel * t)
        dv = net.dist[nbB]
        if 0 <= dv <= t:
            hit = nbB; break
        t += 1
    net.reset(m)
    if hit < 0:
        return np.nan
    m = net.run_bfs(hit, tmax)
    t2 = t + 1; res = np.nan
    while t2 <= t + tmax:
        nA = net.node(pA + vel * t2)
        dv = net.dist[nA]
        if 0 <= dv <= t2 - t:
            res = t2; break
        t2 += 1
    net.reset(m)
    return res


def perp(u, rng):
    a = rng.normal(size=3); a -= (a @ u) * u
    return a / np.linalg.norm(a)


def clock_scan(net, u, c, L, fs, nsamp, rng, nfix=None):
    starts = []
    for k in range(nsamp):
        nvec = perp(u, rng) if nfix is None else nfix
        p = rng.uniform(0, net.Lb, 3)
        starts.append((p, p + L * nvec))
    raw = {}
    for f in fs:
        vel = f * c * u
        tmax = int(3.0 * L / (c * np.sqrt(1 - f * f))) + 30
        raw[f] = [float(tick(net, a, b, vel, tmax)) for a, b in starts]
    bs = np.random.default_rng(5)
    a0 = np.array(raw[fs[0]])
    rate, err = {}, {}
    for f in fs:
        a1 = np.array(raw[f])
        rate[f] = float(np.nanmean(a0) / np.nanmean(a1))
        idx = [bs.integers(0, len(a0), len(a0)) for _ in range(400)]
        err[f] = float(np.std([np.nanmean(a0[i]) / np.nanmean(a1[i]) for i in idx]))
    return rate, err, raw


def tick_blob(net, pA, pB, vel, tmax, rho):
    """Mirrors are small patterns: all nodes within rho of the mirror point (same at rest and moving)."""
    a0 = net.node(pA)
    m = net.run_bfs(a0, tmax)
    t = 1; hit = -1
    while t <= tmax:
        ids = net.tree.query_ball_point((pB + vel * t) % net.Lb, rho)
        if ids:
            dv = net.dist[ids]; ok = (dv >= 0) & (dv <= t)
            if ok.any():
                hit = ids[int(np.argmin(np.where(ok, dv, 10 ** 9)))]; break
        t += 1
    net.reset(m)
    if hit < 0:
        return np.nan
    m = net.run_bfs(hit, tmax)
    t2 = t + 1; res = np.nan
    while t2 <= t + tmax:
        ids = net.tree.query_ball_point((pA + vel * t2) % net.Lb, rho)
        if ids:
            dv = net.dist[ids]
            if ((dv >= 0) & (dv <= t2 - t)).any():
                res = t2; break
        t2 += 1
    net.reset(m)
    return res
