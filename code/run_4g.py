"""Experiment 4g: 0.8c particles on 2x2x2-tiled raw glass s3 (13.8M nodes, 652 A box) and tiled diamond.
Same rule/measures as 4f; lean memory + compiled tracking."""
import numpy as np, json, time, sys, gc
from numba import njit, prange
import scipy.sparse as sp
from scipy.sparse.linalg import cg, LinearOperator
from exp4 import load_net, omega, mu_of_k, step, bipartite_sign
import exp2

eta = 0.15; T = 600; TS = 250; TH = 425
CT = {'glass': 1.5495, 'diamond': 1.3562, 'raw_tiled_test': 1.529}
D2 = {'100': np.array([1., 0, 0]), '123': np.array([1., 2, 3]) / np.sqrt(14)}
STARTS = [(0.5, 0.5, 0.5), (0.25, 0.25, 0.75), (0.75, 0.25, 0.25), (0.25, 0.75, 0.25), (0.75, 0.75, 0.75), (0.5, 0.25, 0.5)]

def load_tiled(name):
    """Tile a base network 2x2x2. rev and bond vectors are copies of the base (tiling preserves slot structure)."""
    base = {'glass': 'glass', 'diamond': 'diamond', 'raw_tiled_test': 'glass_small_raw'}[name]
    if base == 'glass_small_raw':
        from exp4 import G5; L0 = 30 * 4 * 2.35 / np.sqrt(3)
        p0 = np.load(G5 + 'glass_n216000_s1_raw_pos.npy') % L0; n0 = np.load(G5 + 'glass_n216000_s1_raw_nb.npy').astype(np.int64)
        N0 = len(n0); r0 = np.argmax(n0[n0] == np.arange(N0)[:, None, None], axis=2).astype(np.int64)
        d0 = p0[n0] - p0[:, None, :]; d0 -= L0 * np.round(d0 / L0); LB0 = L0
    else:
        p0, n0, r0, d0, LB0 = load_net(base)
    bip0 = bipartite_sign(n0) if name == 'diamond' else None
    P, NB, LB = exp2.tile(p0 % LB0, n0, LB0, 2)
    N0 = len(n0); NB = NB.astype(np.int32); REV = np.tile(r0, (8, 1)).astype(np.int32)
    ok = True
    for s in range(0, len(NB), 2_000_000):
        e = min(s + 2_000_000, len(NB)); ok &= bool((NB[NB[s:e], REV[s:e]] == np.arange(s, e)[:, None]).all())
    assert ok, 'tiled reverse-slot table inconsistent'
    bip = None
    if bip0 is not None:
        bip = np.tile(bip0, 8)
        assert (bip[NB] == -bip[:, None]).all(), 'tiled bipartite signs inconsistent'
    del n0, r0; gc.collect()
    return P, NB, REV, d0, LB, LB0, bip

def adj(nb):
    N = len(nb)
    return sp.csr_matrix((np.ones(4 * N, np.float64), nb.ravel().astype(np.int32), np.arange(0, 4 * N + 1, 4, dtype=np.int64)), shape=(N, N))

@njit(parallel=True, cache=True)
def link_from_xy(x, y, nb, out):
    for i in prange(nb.shape[0]):
        for a in range(4):
            out[i, a] = x[i] + y[nb[i, a]]

@njit(cache=True)
def in_err(phi, nb, X):
    N = nb.shape[0]; acc = np.zeros(N, np.complex128)
    for i in range(N):
        for a in range(4):
            acc[nb[i, a]] += phi[i, a]
    m = 0.0
    for i in range(N):
        v = abs(acc[i] - X[i])
        if v > m: m = v
    return m

def build_packet(pos, nb, d0, LB, A, c0, kvec, vvec, w, sigma, bip=None):
    """Same construction as exp4.build_packet, written to free memory as it goes."""
    rr0 = pos - c0; rr0 -= LB * np.round(rr0 / LB)
    def S_at(t):
        r2 = ((rr0 - vvec * t) ** 2).sum(1)
        return np.exp(-r2 / (4 * sigma ** 2)) * np.exp(1j * (rr0 @ kvec + w * t))
    S0 = S_at(0.0); del rr0; gc.collect()
    ce = np.cos(eta); nC = 400
    xs = np.cos(np.pi * (np.arange(nC) + 0.5) / nC); fv = np.sqrt(1 - ce ** 2 * xs ** 2)
    coef = np.array([2.0 / nC * (fv * np.cos(j * np.pi * (np.arange(nC) + 0.5) / nC)).sum() for j in range(nC)]); coef[0] /= 2
    nuse = int(np.nonzero(np.abs(coef) > 1e-13)[0].max()) + 1
    B = lambda v: (A @ v) * 0.25
    T0, T1 = S0, B(S0); acc = coef[0] * T0 + coef[1] * T1
    for j in range(2, nuse):
        T2 = B(T1); T2 *= 2; T2 -= T0; T0, T1 = T1, T2; acc += coef[j] * T1
    del T0, T1; gc.collect()
    AS0 = A @ S0
    S1 = ce * 0.25 * AS0 + 1j * acc; del acc
    X = np.exp(1j * eta) * (np.cos(eta) / 2 * AS0 - S1); del AS0, S1; gc.collect()
    rhs = 4 * S0 - A @ X
    if bip is not None:
        rhs = rhs - bip * (bip @ rhs) / len(bip)
    M = LinearOperator(A.shape, matvec=lambda v: 16 * v - A @ (A @ v), dtype=complex)
    lam = 4 - mu_of_k(d0, kvec); Xc = np.exp(1j * eta) * (np.cos(eta) / 2 * lam - np.exp(1j * w)); den = 16 - lam ** 2
    p = 1 / 8 if den < 1e-9 else (4 - lam * Xc) / den
    x, info = cg(M, rhs, x0=p * S0, rtol=1e-9, maxiter=5000); del rhs; gc.collect()
    y = (X - A @ x) / 4
    phi = np.empty((len(nb), 4), np.complex128); link_from_xy(x, y, nb, phi); del x, y; gc.collect()
    e_out = float(np.abs(phi.sum(1) - S0).max()); e_in = float(in_err(phi, nb, X))
    return phi, e_in, e_out, info

@njit(parallel=True, cache=True)
def moments(phi, pos, cen, LB, rcut2, kv, off):
    """Weighted nearest-image moments about cen. Returns w, wr(3), w_in, wr_in(3), sum rho r^2, phase sum (re, im)."""
    N = phi.shape[0]
    w = 0.0; wx = 0.0; wy = 0.0; wz = 0.0; v = 0.0; vx = 0.0; vy = 0.0; vz = 0.0; q = 0.0; sre = 0.0; sim = 0.0
    for i in prange(N):
        rho = 0.0; S = 0j
        for a in range(4):
            z = phi[i, a]; rho += z.real * z.real + z.imag * z.imag; S += z
        rx = pos[i, 0] - cen[0]; rx -= LB * np.round(rx / LB)
        ry = pos[i, 1] - cen[1]; ry -= LB * np.round(ry / LB)
        rz = pos[i, 2] - cen[2]; rz -= LB * np.round(rz / LB)
        r2 = rx * rx + ry * ry + rz * rz
        w += rho; wx += rho * rx; wy += rho * ry; wz += rho * rz; q += rho * r2
        if r2 < rcut2:
            v += rho; vx += rho * rx; vy += rho * ry; vz += rho * rz
        ph = -(kv[0] * (rx + off[0]) + kv[1] * (ry + off[1]) + kv[2] * (rz + off[2]))
        e = np.cos(ph) + 1j * np.sin(ph); t = S * e
        sre += t.real; sim += t.imag
    return np.array([w, wx, wy, wz, v, vx, vy, vz, q, sre, sim])

def one(pos, nb, rev, d0, LB, A, bip, c, vf, dn, SIG, c0, Tn=T):
    u = D2[dn]; g = 1 / np.sqrt(1 - vf ** 2); kv = g * eta * vf / c * u; kn = np.linalg.norm(kv)
    w = omega(d0, kv, eta); dk = 1e-4 * u; vg = (omega(d0, kv + dk, eta) - omega(d0, kv - dk, eta)) / 2e-4
    t0 = time.time()
    phi, e_in, e_out, info = build_packet(pos, nb, d0, LB, A, c0, kv, -vg * u, w, SIG, bip)
    tb = time.time() - t0
    out = np.zeros_like(phi); ce = np.cos(eta) / 2; em = np.exp(-1j * eta); z3 = np.zeros(3); rc2 = (0.45 * LB) ** 2
    cen = c0.copy(); un = c0.copy(); cw = c0.copy(); uw = c0.copy(); ts = []; xs = []; xw = []; ph = []; nm = []; rms = []
    for t in range(Tn + 1):
        if t % 5 == 0:
            for _ in range(3):
                m = moments(phi, pos, cen, LB, rc2, kv, z3); sh = m[1:4] / m[0]; cen = (cen + sh) % LB; un = un + sh
            for _ in range(3):
                m = moments(phi, pos, cw, LB, rc2, kv, z3); sh = m[5:8] / m[4]; cw = (cw + sh) % LB; uw = uw + sh
            xw.append(float((uw - c0) @ (-u)))
            m = moments(phi, pos, cen, LB, rc2, kv, un - c0)
            ph.append(float(np.angle(m[9] + 1j * m[10]))); ts.append(t); xs.append(float((un - c0) @ (-u))); nm.append(float(m[0]))
            rms.append(float(np.sqrt(m[8] / m[0])))
        if t < Tn: step(phi, nb, rev, ce, em, out); phi, out = out, phi
    del phi, out; gc.collect()
    return dict(ts=ts, xs=xs, xw=xw, ph=ph, nm=nm, rms=rms, kn=float(kn), e_in=e_in, e_out=e_out, cg_info=int(info), t_build=tb, t_total=time.time() - t0)

def score(r, c, vf, dn, SIG, s0):
    t = np.array(r['ts'], float); x = np.array(r['xs']); m = t >= TS
    v = np.polyfit(t[m], x[m], 1)[0]
    h1 = (t >= TS) & (t <= TH); h2 = t >= TH
    vw = np.polyfit(t[m], np.array(r['xw'])[m], 1)[0]; v1 = np.polyfit(t[h1], x[h1], 1)[0]; v2 = np.polyfit(t[h2], x[h2], 1)[0]
    wc = []; wv = []
    for a0 in range(0, T - 39, 20):
        mm = (t >= a0) & (t <= a0 + 40); wv.append(np.polyfit(t[mm], x[mm], 1)[0]); wc.append(a0 + 20)
    wc = np.array(wc, float); wv = np.array(wv); ms = wc - 20 >= TS
    p = np.polyfit(wc[ms], wv[ms], 1); D = 100 * p[0] * 100 / v; W = 100 * np.sqrt(((wv[ms] - np.polyval(p, wc[ms])) ** 2).mean()) / v
    st = [wc[i] for i in range(len(wv)) if abs(wv[i] / v - 1) < 0.005]; tset = float(st[0]) if st else None
    wm = np.polyfit(t[m], np.unwrap(r['ph'])[m], 1)[0]; clock = wm - r['kn'] * v
    err = 100 * (clock / eta / np.sqrt(1 - (v / c) ** 2) - 1)
    return dict(v_win=float(vw / c), track_diff=float(100 * (vw / v - 1)), vf=vf, dir=dn, sig=SIG, start=s0, v=float(v / c), v1=float(v1 / c), v2=float(v2 / c),
                half=float(100 * (v2 / v1 - 1)), D=float(D), W=float(W), t_settle=tset, win_t=wc.tolist(), win_v=(wv / c).tolist(), clock_err=float(err),
                norm=float(abs(r['nm'][-1] / r['nm'][0] - 1)), rms_end=r['rms'][-1], e_in=r['e_in'], e_out=r['e_out'], cg_info=r['cg_info'],
                t_build=r['t_build'], t_total=r['t_total'])

if __name__ == '__main__':
    res = []
    for name in ['diamond', 'glass']:
        t0 = time.time(); pos, nb, rev, d0, LB, LB0, bip = load_tiled(name); A = adj(nb); c = CT[name]
        print(f'{name}: N={len(nb)} LB={LB:.1f} loaded {time.time()-t0:.0f}s', flush=True)
        jobs = [('100', STARTS[0]), ('123', STARTS[0])] if name == 'diamond' else [(dn, s) for dn in ['100', '123'] for s in STARTS]
        for dn, s in jobs:
            c0 = np.array(s) * LB0          # same absolute positions (A) as Experiment 4f
            r = one(pos, nb, rev, d0, LB, A, bip, c, 0.8, dn, 40.0, c0)
            sc = score(r, c, 0.8, dn, 40.0, s); sc['net'] = name; res.append(sc)
            json.dump(res, open('res_4g.json', 'w'))
            print(f"{name} {dn} start{s}: v/c {sc['v']:.4f} halves {sc['half']:+.2f}% D {sc['D']:+.2f}%/100t W {sc['W']:.2f}% trackdiff {sc['track_diff']:+.2f}% "
                  f"clock {sc['clock_err']:+.2f}% norm {sc['norm']:.1e} rms_end {sc['rms_end']:.0f} e_out {sc['e_out']:.1e} cg {sc['cg_info']} "
                  f"build {sc['t_build']:.0f}s total {sc['t_total']:.0f}s", flush=True)
        del pos, nb, rev, A; gc.collect()
    print('DONE', flush=True)
