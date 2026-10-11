"""Experiment 4: GBU-native amplitude rule (selection with memory, direction-valued, phase-carrying).
State phi[i,a]: amplitude at node i heading out along link a.  Tick: mix (coin) then move (flip-flop shift).
"""
import numpy as np, json, sys, time
from numba import njit, prange
import scipy.sparse as sp
from scipy.sparse.linalg import cg, LinearOperator
import os
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
G5 = os.environ.get('GBU_NETWORKS', os.path.join(HERE, '..', 'networks')) + '/'

def load_net(name):
    if name == 'glass':
        pos = np.load(G5 + 'glass_n1728000_s3_raw_pos.npy'); nb = np.load(G5 + 'glass_n1728000_s3_raw_nb.npy').astype(np.int64)
        LB = 2 * 30 * 4 * 2.35 / np.sqrt(3)
    elif name == 'glass4':   # independent replication sample (Exp 4i-R), same recipe as s3, seed 4
        pos = np.load(G5 + 'glass_n1728000_s4_raw_pos.npy'); nb = np.load(G5 + 'glass_n1728000_s4_raw_nb.npy').astype(np.int64)
        LB = 2 * 30 * 4 * 2.35 / np.sqrt(3)
    elif name in ('smooth_tiled', 'raw_tiled'):
        import exp2
        tag = 's1_final' if name == 'smooth_tiled' else 's1_raw'
        L0 = 30 * 4 * 2.35 / np.sqrt(3)
        p0 = np.load(G5 + f'glass_n216000_{tag}_pos.npy'); n0 = np.load(G5 + f'glass_n216000_{tag}_nb.npy').astype(np.int64)
        pos, nb, LB = exp2.tile(p0 % L0, n0, L0, 2)
    elif name == 'glass_small':
        pos = np.load(G5 + 'glass_n216000_s1_final_pos.npy'); nb = np.load(G5 + 'glass_n216000_s1_final_nb.npy').astype(np.int64)
        LB = 30 * 4 * 2.35 / np.sqrt(3)
    else:
        import crn
        ncell = 60 if name == 'diamond' else 30
        pos, nb, LB = crn.diamond(ncell); nb = nb.astype(np.int64)
    pos = pos % LB
    N = len(nb)
    # reverse slot: rev[i,a] = b such that nb[nb[i,a], b] == i
    j = nb; rev = np.argmax(nb[j] == np.arange(N)[:, None, None], axis=2).astype(np.int64)
    assert (nb[j, rev] == np.arange(N)[:, None]).all()
    d = pos[nb] - pos[:, None, :]; d -= LB * np.round(d / LB)
    return pos, nb, rev, d, LB

def mu_of_k(d, kvec):
    """Effective Laplacian eigenvalue for a plane wave: mean over nodes of sum_a (1 - cos k.d)."""
    return float((1 - np.cos(d @ kvec)).sum(1).mean())

def omega(d, kvec, eta):
    return float(np.arccos(np.cos(eta) * (1 - mu_of_k(d, kvec) / 4)))

@njit(parallel=True, cache=True)
def step(phi, nb, rev, ce, em, out):
    N = phi.shape[0]
    for i in prange(N):
        s = phi[i, 0] + phi[i, 1] + phi[i, 2] + phi[i, 3]
        for a in range(4):
            out[nb[i, a], rev[i, a]] = ce * s - em * phi[i, a]

@njit(parallel=True, cache=True)
def step_classical(p, nb, rev, out):
    # squared weights of the Grover-type mixing at eta=0: keep-back 1/4, each other 1/4 -> uniform
    N = p.shape[0]
    for i in prange(N):
        s = 0.25 * (p[i, 0] + p[i, 1] + p[i, 2] + p[i, 3])
        for a in range(4):
            out[nb[i, a], rev[i, a]] = s

def adj(nb):
    N = len(nb); r = np.repeat(np.arange(N), 4)
    return sp.csr_matrix((np.ones(4 * N), (r, nb.ravel())), shape=(N, N))

def build_packet(pos, nb, d, LB, A, c0, kvec, vvec, w, eta, sigma, bip=None):
    """Link amplitudes whose node totals are S0 now and S1 next tick (branch S ~ e^{+i w t}, moving along vvec).
    Least-norm solution phi(i,a) = x_i + y_{nb[i,a]}; solved by CG starting from the local plane-wave ansatz."""
    rr0 = pos - c0; rr0 -= LB * np.round(rr0 / LB)
    def S_at(t):
        r = rr0 - vvec * t
        return np.exp(-(r ** 2).sum(1) / (4 * sigma ** 2)) * np.exp(1j * (rr0 @ kvec + w * t))
    S0 = S_at(0.0)
    # exact one-branch next state: S1 = e^{+i Omega(A)} S0,  cos Omega = cos(eta) A/4
    #   = (cos eta/4) A S0 + i sqrt(I - cos^2(eta) A^2/16) S0   (Chebyshev series for the square root)
    ce = np.cos(eta); nC = 400
    xs = np.cos(np.pi * (np.arange(nC) + 0.5) / nC); fv = np.sqrt(1 - ce ** 2 * xs ** 2)
    coef = np.array([2.0 / nC * (fv * np.cos(j * np.pi * (np.arange(nC) + 0.5) / nC)).sum() for j in range(nC)]); coef[0] /= 2
    nuse = int(np.nonzero(np.abs(coef) > 1e-13)[0].max()) + 1
    B = lambda v: (A @ v) / 4
    T0, T1 = S0, B(S0); acc = coef[0] * T0 + coef[1] * T1
    for j in range(2, nuse):
        T0, T1 = T1, 2 * B(T1) - T0; acc = acc + coef[j] * T1
    S1 = ce * B(S0) + 1j * acc
    X = np.exp(1j * eta) * (np.cos(eta) / 2 * (A @ S0) - S1)
    rhs = 4 * S0 - A @ X
    if bip is not None:                      # remove the (unused) staggered mode of a bipartite crystal
        rhs = rhs - bip * (bip @ rhs) / len(bip)
    M = LinearOperator(A.shape, matvec=lambda v: 16 * v - A @ (A @ v), dtype=complex)
    lam = 4 - mu_of_k(d, kvec); Xc = np.exp(1j * eta) * (np.cos(eta) / 2 * lam - np.exp(1j * w)); den = 16 - lam ** 2
    p = 1 / 8 if den < 1e-9 else (4 - lam * Xc) / den
    x, info = cg(M, rhs, x0=p * S0, rtol=1e-9, maxiter=5000)
    y = (X - A @ x) / 4
    phi = x[:, None] + y[nb]
    e_out = np.abs(phi.sum(1) - S0).max()
    e_in = np.abs(np.bincount(nb.ravel(), weights=phi.real.ravel(), minlength=len(nb)) + 1j * np.bincount(nb.ravel(), weights=phi.imag.ravel(), minlength=len(nb)) - X).max()
    return phi.astype(np.complex128), float(e_in), float(e_out)

def run(name, pos, nb, rev, d, LB, A, eta, kvec, sigma, T, c_est, every=5, classical=False, bip=None, filt=0, dens='link'):
    c0 = np.array([LB / 2] * 3)
    w = omega(d, kvec, eta)
    kn = np.linalg.norm(kvec)
    if kn > 0:
        dk = 1e-4 * kvec / kn
        vg = (omega(d, kvec + dk, eta) - omega(d, kvec - dk, eta)) / (2e-4)
        vvec = -vg * kvec / kn          # + branch moves against k
    else:
        vvec = np.zeros(3)
    phi, e_in, e_out = build_packet(pos, nb, d, LB, A, c0, kvec, vvec, w, eta, sigma, bip)
    out = np.zeros_like(phi); ce = np.cos(eta) / 2; em = np.exp(-1j * eta)
    hid = 0.0
    if filt:   # remove the motionless hidden modes (link patterns with zero node totals): U^2 = e^{-2i eta} on them
        n0 = np.vdot(phi, phi).real
        a = phi.copy(); step(a, nb, rev, ce, em, out); b = np.zeros_like(a); step(out.copy(), nb, rev, ce, em, b)
        for _ in range(filt):
            a = phi.copy(); o1 = np.zeros_like(a); step(a, nb, rev, ce, em, o1); o2 = np.zeros_like(a); step(o1, nb, rev, ce, em, o2)
            phi = (o2 - np.exp(-2j * eta) * phi) / (np.exp(2j * w) - np.exp(-2j * eta))
        hid = float(1 - np.vdot(phi, phi).real / n0)
    if classical:
        p = np.abs(phi) ** 2; outp = np.zeros_like(p)
    rec = dict(t=[], c=[], ph=[], phw=[], rms=[], fout=[], norm=[])
    cen = c0.copy(); un = c0.copy()
    for t in range(T + 1):
        if t % every == 0:
            rho = ((np.abs(phi) ** 2).sum(1) if dens == 'link' else np.abs(phi.sum(1)) ** 2) if not classical else p.sum(1)
            for _ in range(2):
                r = pos - cen; r -= LB * np.round(r / LB)
                sel = (r ** 2).sum(1) < (3.5 * sigma + 0.5 * t) ** 2
                shift = (rho[sel, None] * r[sel]).sum(0) / rho[sel].sum()
                cen = (cen + shift) % LB; un = un + shift
            r = pos - cen; r -= LB * np.round(r / LB); r2 = (r ** 2).sum(1)
            rms = np.sqrt((rho * r2).sum() / rho.sum())
            fout = float(rho[r2 > (3 * rms) ** 2].sum() / rho.sum())
            if not classical:
                S = phi.sum(1); R = r + (un - c0)          # continuous coordinate from the start point
                rec['ph'].append(float(np.angle((S * np.exp(-1j * (R @ kvec))).sum())))   # carrier (Fourier) mode
                win = np.exp(-r2 / (2 * sigma ** 2))
                rec['phw'].append(float(np.angle((win * S * np.exp(-1j * (r @ kvec))).sum())))  # centre phase
            rec['t'].append(t); rec['c'].append(un.tolist()); rec['rms'].append(float(rms)); rec['fout'].append(fout)
            rec['norm'].append(float(rho.sum()))
        if t < T:
            if classical:
                step_classical(p, nb, rev, outp); p, outp = outp, p
            else:
                step(phi, nb, rev, ce, em, out); phi, out = out, phi
    return dict(rec=rec, w=w, vg=float(np.linalg.norm(vvec)), e_in=e_in, e_out=e_out, hid=hid)

def analyse(res, kvec):
    rec = res['rec']; t = np.array(rec['t'], float); c = np.array(rec['c'])
    kn = np.linalg.norm(kvec); n = -kvec / kn if kn > 0 else np.array([1.0, 0, 0])
    x = (c - c[0]) @ n
    m = t >= 20; h = m & (t < t[m].mean()); h2 = m & (t >= t[m].mean())
    v = np.polyfit(t[m], x[m], 1)[0]; v1 = np.polyfit(t[h], x[h], 1)[0]; v2 = np.polyfit(t[h2], x[h2], 1)[0]
    out = dict(v=float(v), v1=float(v1), v2=float(v2), norm_drift=float(abs(rec['norm'][-1] / rec['norm'][0] - 1)),
               fout_end=rec['fout'][-1], rms0=rec['rms'][0], rms_end=rec['rms'][-1])
    if rec['ph']:
        ph = np.unwrap(np.array(rec['ph']))
        wm = float(np.polyfit(t[m], ph[m], 1)[0])          # measured carrier frequency
        out['w_meas'] = wm; out['clock'] = float(wm - kn * v)   # phase rate at the moving particle
        phw = np.unwrap(np.array(rec['phw'])); out['clock_centre'] = float(abs(np.polyfit(t[m], phw[m], 1)[0]))
    return out

def bipartite_sign(nb):
    """Return +-1 sublattice signs if the network is bipartite, else None."""
    N = len(nb); col = np.zeros(N, int); col[0] = 1; stack = [0]
    while stack:
        i = stack.pop()
        for j in nb[i]:
            if col[j] == 0: col[j] = -col[i]; stack.append(j)
            elif col[j] == col[i]: return None
    return col.astype(float)
