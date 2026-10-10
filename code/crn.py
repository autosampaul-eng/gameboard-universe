"""Continuous random network (glass) builder: WWW bond switching with the Keating
potential, local + global relaxation, Metropolis acceptance (Barkema-Mousseau style).
Every node keeps exactly 4 links. Periodic cubic box. Units: Angstrom, eV."""
import numpy as np
from numba import njit

D = 2.35          # ideal bond length
ALPHA = 2.965     # bond-stretch constant
BETA = 0.285 * ALPHA
C1 = 3.0 * ALPHA / (16.0 * D * D)
C2 = 3.0 * BETA / (8.0 * D * D)
D2_3 = D * D / 3.0


def diamond(ncell):
    a = 4.0 * D / np.sqrt(3.0)
    basis = np.array([[0, 0, 0], [0, .5, .5], [.5, 0, .5], [.5, .5, 0],
                      [.25, .25, .25], [.25, .75, .75], [.75, .25, .75], [.75, .75, .25]])
    g = np.array([[i, j, k] for i in range(ncell) for j in range(ncell) for k in range(ncell)], float)
    pos = (g[:, None, :] + basis[None, :, :]).reshape(-1, 3) * a
    Lb = ncell * a
    from scipy.spatial import cKDTree
    t = cKDTree(pos, boxsize=Lb)
    pairs = t.query_pairs(D * 1.05, output_type='ndarray')
    n = len(pos)
    nb = -np.ones((n, 4), np.int64); cnt = np.zeros(n, np.int64)
    for i, j in pairs:
        nb[i, cnt[i]] = j; cnt[i] += 1; nb[j, cnt[j]] = i; cnt[j] += 1
    assert (cnt == 4).all()
    return pos, nb, Lb


@njit(cache=True)
def mi(x, Lb):
    return x - Lb * np.round(x / Lb)


@njit(cache=True)
def atom_e(i, pos, nb, Lb):
    e = 0.0
    r = np.empty((4, 3))
    for a in range(4):
        j = nb[i, a]
        for c in range(3):
            r[a, c] = mi(pos[j, c] - pos[i, c], Lb)
        s = r[a, 0] ** 2 + r[a, 1] ** 2 + r[a, 2] ** 2 - D * D
        e += 0.5 * C1 * s * s
    for a in range(4):
        for b in range(a + 1, 4):
            p = r[a, 0] * r[b, 0] + r[a, 1] * r[b, 1] + r[a, 2] * r[b, 2] + D2_3
            e += C2 * p * p
    return e


@njit(cache=True)
def atom_eg(i, pos, nb, Lb, g):
    e = 0.0
    r = np.empty((4, 3))
    for a in range(4):
        j = nb[i, a]
        for c in range(3):
            r[a, c] = mi(pos[j, c] - pos[i, c], Lb)
        s = r[a, 0] ** 2 + r[a, 1] ** 2 + r[a, 2] ** 2 - D * D
        e += 0.5 * C1 * s * s
        for c in range(3):
            f = 0.5 * C1 * 4.0 * s * r[a, c]
            g[j, c] += f; g[i, c] -= f
    for a in range(4):
        ja = nb[i, a]
        for b in range(a + 1, 4):
            jb = nb[i, b]
            p = r[a, 0] * r[b, 0] + r[a, 1] * r[b, 1] + r[a, 2] * r[b, 2] + D2_3
            e += C2 * p * p
            for c in range(3):
                fa = 2.0 * C2 * p * r[b, c]; fb = 2.0 * C2 * p * r[a, c]
                g[ja, c] += fa; g[jb, c] += fb; g[i, c] -= fa + fb
    return e


@njit(cache=True)
def total_e(pos, nb, Lb):
    e = 0.0
    for i in range(pos.shape[0]):
        e += atom_e(i, pos, nb, Lb)
    return e


@njit(cache=True)
def gather(seeds, nseed, radius, nb, stamp, tag, out, maxout):
    """BFS up to `radius` from seeds; returns count. stamp[] marks membership with tag."""
    n = 0
    for s in range(nseed):
        x = seeds[s]
        if stamp[x] != tag:
            stamp[x] = tag; out[n] = x; n += 1
    lo = 0
    for _ in range(radius):
        hi = n
        for q in range(lo, hi):
            x = out[q]
            for a in range(4):
                y = nb[x, a]
                if stamp[y] != tag and n < maxout:
                    stamp[y] = tag; out[n] = y; n += 1
        lo = hi
    return n


@njit(cache=True)
def fire(pos, nb, Lb, mov, nmov, cen, ncen, g, ftol, maxit, Ecut=1e300):
    """FIRE minimisation of atoms mov[:nmov]; energy from centres cen[:ncen]."""
    dt = 0.05; dtmax = 0.3; alpha = 0.1; npos = 0
    v = np.zeros((nmov, 3))
    E = 0.0
    for it in range(maxit):
        for q in range(ncen):
            x = cen[q]
            for c in range(3):
                g[x, c] = 0.0
            for a in range(4):
                y = nb[x, a]
                for c in range(3):
                    g[y, c] = 0.0
        E = 0.0
        for q in range(ncen):
            E += atom_eg(cen[q], pos, nb, Lb, g)
        P = 0.0; vn = 0.0; fn = 0.0; fmax = 0.0
        for q in range(nmov):
            x = mov[q]
            for c in range(3):
                f = -g[x, c]
                P += f * v[q, c]; vn += v[q, c] ** 2; fn += f * f
                if abs(f) > fmax:
                    fmax = abs(f)
        if fmax < ftol:
            break
        if it == 60 and E > Ecut + 1.0:
            break
        if it == 150 and E > Ecut + 0.3:
            break
        vn = np.sqrt(vn); fn = np.sqrt(fn) + 1e-30
        if P > 0:
            for q in range(nmov):
                x = mov[q]
                for c in range(3):
                    v[q, c] = (1 - alpha) * v[q, c] + alpha * vn * (-g[x, c]) / fn
            npos += 1
            if npos > 5:
                dt = min(dt * 1.1, dtmax); alpha *= 0.99
        else:
            npos = 0; dt *= 0.5; alpha = 0.1
            for q in range(nmov):
                for c in range(3):
                    v[q, c] = 0.0
        for q in range(nmov):
            x = mov[q]
            for c in range(3):
                v[q, c] += -g[x, c] * dt
                pos[x, c] += v[q, c] * dt
    return E


@njit(cache=True)
def has_edge(nb, x, y):
    for a in range(4):
        if nb[x, a] == y:
            return True
    return False


@njit(cache=True)
def short_ring(nb, x, y, stamp, tag, buf):
    """True if y is reachable from x within 3 steps without the direct x-y edge
    (i.e. the x-y bond closes a ring of size <= 4)."""
    stamp[x] = tag; buf[0] = x; n = 1; lo = 0
    for d in range(3):
        hi = n
        for q in range(lo, hi):
            u = buf[q]
            for a in range(4):
                w = nb[u, a]
                if (u == x and w == y) or (u == y and w == x):
                    continue
                if w == y:
                    return True
                if stamp[w] != tag:
                    stamp[w] = tag; buf[n] = w; n += 1
        lo = hi
    return False


@njit(cache=True)
def replace(nb, x, old, new):
    for a in range(4):
        if nb[x, a] == old:
            nb[x, a] = new
            return


@njit(cache=True)
def mc(pos, nb, Lb, nsteps, kT, radius, seed, stats):
    np.random.seed(seed)
    n = pos.shape[0]
    stamp = np.zeros(n, np.int64); tag = 1
    mov = np.empty(4096, np.int64); cen = np.empty(8192, np.int64)
    buf = np.empty(512, np.int64)
    g = np.zeros((n, 3))
    save = np.empty((4096, 3))
    seeds = np.empty(4, np.int64)
    acc = 0; tried = 0
    for step in range(nsteps):
        B = np.random.randint(n)
        ib = np.random.randint(4)
        C = nb[B, ib]
        A = nb[B, (ib + 1 + np.random.randint(3)) % 4]
        Dd = -1
        o = np.random.randint(4)
        for t in range(4):
            cand_d = nb[C, (o + t) % 4]
            if cand_d != B and cand_d != A and not has_edge(nb, cand_d, B):
                Dd = cand_d; break
        if Dd < 0 or has_edge(nb, A, C):
            continue
        tried += 1
        # region (old topology already contains A,B,C,D; union with new is the same nodes)
        seeds[0] = A; seeds[1] = B; seeds[2] = C; seeds[3] = Dd
        tag += 1
        nmov = gather(seeds, 4, radius, nb, stamp, tag, mov, 4096)
        tag += 1
        ncen = gather(mov, nmov, 1, nb, stamp, tag, cen, 8192)
        Eold = 0.0
        for q in range(ncen):
            Eold += atom_e(cen[q], pos, nb, Lb)
        for q in range(nmov):
            for c in range(3):
                save[q, c] = pos[mov[q], c]
        # switch: A-B, C-D  ->  A-C, B-D
        replace(nb, A, B, C); replace(nb, B, A, Dd); replace(nb, C, Dd, A); replace(nb, Dd, C, B)
        tag += 1
        bad = short_ring(nb, A, C, stamp, tag, buf)
        if not bad:
            tag += 1
            bad = short_ring(nb, B, Dd, stamp, tag, buf)
        ok = False
        if not bad:
            thr = -kT * np.log(np.random.random() + 1e-300)
            Enew = fire(pos, nb, Lb, mov, nmov, cen, ncen, g, 2e-3, 400, Eold + thr)
            if Enew - Eold < thr:
                ok = True
        if ok:
            acc += 1
        else:
            replace(nb, A, C, B); replace(nb, B, Dd, A); replace(nb, C, A, Dd); replace(nb, Dd, B, C)
            for q in range(nmov):
                for c in range(3):
                    pos[mov[q], c] = save[q, c]
    stats[0] = tried; stats[1] = acc


@njit(cache=True)
def global_relax(pos, nb, Lb, ftol, maxit):
    n = pos.shape[0]
    mov = np.arange(n); g = np.zeros((n, 3))
    return fire(pos, nb, Lb, mov, n, mov, n, g, ftol, maxit)


def quality(pos, nb, Lb):
    n = len(pos)
    r = pos[nb] - pos[:, None, :]
    r -= Lb * np.round(r / Lb)
    bl = np.linalg.norm(r, axis=2)
    ang = []
    for a in range(4):
        for b in range(a + 1, 4):
            c = (r[:, a] * r[:, b]).sum(1) / (bl[:, a] * bl[:, b])
            ang.append(np.degrees(np.arccos(np.clip(c, -1, 1))))
    ang = np.concatenate(ang)
    a0 = 4.0 * D / np.sqrt(3.0)
    q = 2 * np.pi / a0 * np.array([1.0, 1.0, 1.0])
    S111 = float(np.abs(np.exp(1j * pos @ q).sum()) ** 2 / n)
    return dict(n=n, bond_rms=float(np.sqrt(((bl - D) ** 2).mean()) / D), angle_mean=float(ang.mean()),
                angle_rms=float(np.sqrt(((ang - 109.4712) ** 2).mean())), S111=S111,
                E_per_atom=float(total_e(pos, nb, Lb) / n))


@njit(cache=True)
def angle_rings(nb, maxd):
    """For every bond angle j-i-k: size of shortest ring through it (d(j,k) avoiding i) + 2.
    Returns histogram of ring sizes (index = size), over all 6N angles."""
    n = nb.shape[0]
    hist = np.zeros(maxd + 3, np.int64)
    dist = -np.ones(n, np.int64); buf = np.empty(n, np.int64)
    for i in range(n):
        for a in range(4):
            j = nb[i, a]
            # BFS from j avoiding i, up to depth maxd
            dist[j] = 0; buf[0] = j; m = 1; lo = 0
            for d in range(maxd):
                hi = m
                for q in range(lo, hi):
                    u = buf[q]
                    for c in range(4):
                        w = nb[u, c]
                        if w == i or dist[w] >= 0:
                            continue
                        dist[w] = d + 1; buf[m] = w; m += 1
                lo = hi
            for b in range(a + 1, 4):
                k = nb[i, b]
                if dist[k] >= 0:
                    hist[dist[k] + 2] += 1
                else:
                    hist[maxd + 2] += 1
            for q in range(m):
                dist[buf[q]] = -1
    return hist


def bragg(pos, Lb):
    a0 = 4.0 * D / np.sqrt(3.0); n = len(pos)
    out = {}
    for name, h in [('111', (1, 1, 1)), ('220', (2, 2, 0)), ('311', (3, 1, 1))]:
        q = 2 * np.pi / a0 * np.array(h, float)
        out[name] = float(np.abs(np.exp(1j * pos @ q).sum()) ** 2 / n)
    return out
