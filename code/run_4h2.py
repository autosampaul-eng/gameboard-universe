"""Experiment 4h (rerun, amendment 1): particle-only tracker weighted by the conserved link probability.
Change vs run_4h.py: before the inverse transform, the forward-cone spectrum is multiplied by sqrt(w(k)) / W_cic(k), where
  w(k)     = (link probability per node) / |node total|^2 for the one-branch plane-wave mode at k (mean over nodes),
  W_cic(k) = the cloud-in-cell window, prod sinc^2(k_i h / 2).
So |psi_f|^2 measures the conserved probability (as the whole-box centre does) instead of |node total|^2.
Scattering measure s(t) unchanged."""
import numpy as np, json, time, sys, os, gc, argparse
from numba import njit, prange
import run_4h as H
from run_4h import cic, grid_centroid, masks, speed_stats, jobs, CT, T, TS, TH
from run_4g import load_tiled, adj, build_packet, moments, D2, STARTS, eta
from exp4 import load_net, omega, step, mu_of_k
import exp2

@njit(parallel=True, cache=True)
def mode_weights(kpts, d, ce, eta):
    """For each k: mean over nodes of link prob / |node total|^2 for the plane-wave ansatz used to build packets."""
    M = kpts.shape[0]; N = d.shape[0]; out = np.empty(M)
    em = np.exp(1j * eta)
    for m in prange(M):
        kx, ky, kz = kpts[m, 0], kpts[m, 1], kpts[m, 2]
        mu = 0.0
        for i in range(N):
            for a in range(4):
                mu += 1 - np.cos(kx * d[i, a, 0] + ky * d[i, a, 1] + kz * d[i, a, 2])
        mu /= N; lam = 4 - mu
        w = np.arccos(ce * (1 - mu / 4))
        Xc = em * (ce / 2 * lam - np.exp(1j * w)); den = 16 - lam * lam
        p = 1 / 8 if den < 1e-9 else (4 - lam * Xc) / den
        y = (Xc - lam * p) / 4
        num = 0.0; sq = 0.0
        for i in range(N):
            S = 4 * p + 0j
            for a in range(4):
                e = np.cos(kx * d[i, a, 0] + ky * d[i, a, 1] + kz * d[i, a, 2]) + 1j * np.sin(kx * d[i, a, 0] + ky * d[i, a, 1] + kz * d[i, a, 2])
                z = p + y * e; num += z.real * z.real + z.imag * z.imag
                S += y * e
            sq += S.real * S.real + S.imag * S.imag
        out[m] = num / sq
    return out

def tracker_weight(fwd, n, LB, d0, nsub=20000):
    k1 = 2 * np.pi * np.fft.fftfreq(n, d=LB / n)
    idx = np.nonzero(fwd); kp = np.stack([k1[idx[0]], k1[idx[1]], k1[idx[2]]], 1)
    sub = d0[np.random.default_rng(0).choice(len(d0), min(nsub, len(d0)), replace=False)]
    w = mode_weights(kp, np.ascontiguousarray(sub), np.cos(eta), eta)
    h = LB / n; Wc = np.prod(np.sinc(kp * h / 2 / np.pi) ** 2, 1)
    Wt = np.zeros((n, n, n)); Wt[idx] = np.sqrt(w) / Wc
    return Wt, float(w.min()), float(w.max())

def one(pos, nb, rev, d0, LB, A, bip, c, vf, dn, SIG, c0, Tn=T):
    u = D2[dn]; g = 1 / np.sqrt(1 - vf ** 2); kv = g * eta * vf / c * u; kn = np.linalg.norm(kv)
    w = omega(d0, kv, eta); dk = 1e-4 * u; vg = (omega(d0, kv + dk, eta) - omega(d0, kv - dk, eta)) / 2e-4
    t0 = time.time()
    phi, e_in, e_out, info = build_packet(pos, nb, d0, LB, A, c0, kv, -vg * u, w, SIG, bip)
    tb = time.time() - t0
    n = 128 if LB > 500 else (64 if LB > 250 else 32)
    shell, fwd = masks(n, LB, kv); grid = np.zeros((n, n, n), np.complex128)
    Wt, wmin, wmax = tracker_weight(fwd, n, LB, d0)
    out = np.zeros_like(phi); ce = np.cos(eta) / 2; em = np.exp(-1j * eta); z3 = np.zeros(3); rc2 = (0.45 * LB) ** 2
    cen = c0.copy(); un = c0.copy(); cc = c0.copy(); uc = c0.copy(); co = c0.copy(); uo = c0.copy()
    rec = {k: [] for k in ['ts', 'xs', 'xc', 'xo', 'ph', 'nm', 'rms', 'rms_c', 's', 'off', 'fwd_w']}
    for t in range(Tn + 1):
        if t % 5 == 0:
            for _ in range(3):
                m = moments(phi, pos, cen, LB, rc2, kv, z3); sh = m[1:4] / m[0]; cen = (cen + sh) % LB; un = un + sh
            m = moments(phi, pos, cen, LB, rc2, kv, un - c0)
            grid[:] = 0; cic(phi, pos, LB, n, grid); F = np.fft.fftn(grid); P = np.abs(F) ** 2
            Ptot = P.sum(); Psh = P[shell].sum(); Pf = P[fwd].sum()
            rho_c = np.abs(np.fft.ifftn(F * Wt)) ** 2                    # corrected (probability-weighted)
            rho_o = np.abs(np.fft.ifftn(np.where(fwd, F, 0))) ** 2      # original 4h tracker, kept for comparison
            for _ in range(3):
                gm = grid_centroid(rho_c, cc, LB); sh = gm[1:4] / gm[0]; cc = (cc + sh) % LB; uc = uc + sh
                go = grid_centroid(rho_o, co, LB); sh = go[1:4] / go[0]; co = (co + sh) % LB; uo = uo + sh
            gm = grid_centroid(rho_c, cc, LB)
            rec['ts'].append(t); rec['xs'].append(float((un - c0) @ (-u))); rec['xc'].append(float((uc - c0) @ (-u))); rec['xo'].append(float((uo - c0) @ (-u)))
            rec['ph'].append(float(np.angle(m[9] + 1j * m[10]))); rec['nm'].append(float(m[0])); rec['rms'].append(float(np.sqrt(m[8] / m[0])))
            rec['rms_c'].append(float(np.sqrt(gm[4] / gm[0]))); rec['s'].append(float(1 - Pf / Psh)); rec['off'].append(float(1 - Psh / Ptot)); rec['fwd_w'].append(float(Pf / Ptot))
        if t < Tn: step(phi, nb, rev, ce, em, out); phi, out = out, phi
    del phi, out; gc.collect()
    rec.update(kn=float(kn), e_in=e_in, e_out=e_out, cg_info=int(info), t_build=tb, t_total=time.time() - t0, grid_n=n, w_min=wmin, w_max=wmax, vg_carrier=float(vg / c))
    return rec

def score(r, c, Tn=T, ts=TS):
    sc = H.score(r, c, Tn, ts)
    t = np.array(r['ts'], float)
    sc['coh_orig'], _ = speed_stats(t, np.array(r['xo']), c, Tn, ts, (ts + Tn) / 2)
    sc.update(w_min=r['w_min'], w_max=r['w_max'], vg_carrier=r['vg_carrier'])
    return sc

def load(name):
    if name == 'diamond_tiled_test':
        p0, n0, r0, d0, LB0 = load_net('diamond_small')
        from exp4 import bipartite_sign
        bip0 = bipartite_sign(n0)
        P, NB, LB = exp2.tile(p0 % LB0, n0, LB0, 2)
        return P, NB.astype(np.int32), np.tile(r0, (8, 1)).astype(np.int32), d0, LB, LB0, np.tile(bip0, 8)
    return H.load(name)

if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('job', nargs='?', type=int); ap.add_argument('--out', default='../results/exp4h2')
    ap.add_argument('--test'); ap.add_argument('--T', type=int, default=T); ap.add_argument('--sig', type=float, default=40.0); ap.add_argument('--vf', type=float, default=0.8)
    ap.add_argument('--dir', default='100'); a = ap.parse_args()
    spec = dict(net=a.test, vf=a.vf, dir=a.dir, start=STARTS[1], set='test') if a.test else jobs()[a.job]
    t0 = time.time(); pos, nb, rev, d0, LB, LB0, bip = load(spec['net']); A = adj(nb); c = CT.get(spec['net'], 1.3562)
    print(f"job {a.job} {spec}: N={len(nb)} box={LB:.1f} loaded {time.time()-t0:.0f}s", flush=True)
    c0 = np.array(spec['start']) * LB0; Tn = a.T if a.test else T; ts = min(TS, Tn // 2)
    r = one(pos, nb, rev, d0, LB, A, bip, c, spec['vf'], spec['dir'], a.sig if a.test else 40.0, c0, Tn)
    sc = score(r, c, Tn, ts); sc.update(spec); sc['N'] = len(nb); sc['box'] = LB
    print(f"particle-only v/c {sc['coh']['v']:.4f} (orig 4h tracker {sc['coh_orig']['v']:.4f}, whole-box {sc['whole']['v']:.4f}, carrier vg {sc['vg_carrier']:.4f}) "
          f"D {sc['coh']['D']:+.2f} halves {sc['coh']['half']:+.2f}% | ds250 {sc['ds_startup']:+.2f} Gamma {sc['Gamma']:+.3f} max|ds| {sc['ds_max']:.2f} | "
          f"clock {sc['clock_err']:+.2f}% norm {sc['norm']:.1e} w {sc['w_min']:.4f}-{sc['w_max']:.4f} build {sc['t_build']:.0f}s total {sc['t_total']:.0f}s", flush=True)
    if not a.test:
        os.makedirs(a.out, exist_ok=True); json.dump(sc, open(f"{a.out}/job_{a.job:02d}.json", 'w'))
    else:
        json.dump(sc, open(f"test_4h2_{a.test}_{a.dir}_{a.vf}_{a.sig:.0f}.json", 'w'))
