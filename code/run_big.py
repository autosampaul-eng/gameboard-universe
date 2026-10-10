"""Memory-lean runner (same measurements as run_4h2 / run_4i) for tiled boxes, including 3x3x3.
Usage: python run_big.py --net NAME --tile T --point P [--ticks N] [--dir 100] [--out FILE] [--compare]
NAME: glass (1.73M s3 base), diamond (1.73M base), g216 (216k raw glass base), d216 (216k diamond base)."""
import numpy as np, json, time, os, gc, argparse, resource
import lean
from exp4 import load_net, omega, mu_of_k, bipartite_sign
from run_4g import moments, D2, STARTS, eta
from run_4h import cic, grid_centroid, masks, CT, TS
import run_4h2 as H2
from run_4i import POINTS, decay

def rss():
    for l in open('/proc/self/status'):
        if l.startswith('VmRSS'): return int(l.split()[1]) / 1e6
def peak(): return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6

def load(name, T):
    G5 = '/home/user/workspace/gbu5/' if os.path.exists('/home/user/workspace/gbu5') else os.environ.get('GBU_NETWORKS', '../networks') + '/'
    if name == 'g216':
        L0 = 30 * 4 * 2.35 / np.sqrt(3)
        p0 = np.load(G5 + 'glass_n216000_s1_raw_pos.npy') % L0; n0 = np.load(G5 + 'glass_n216000_s1_raw_nb.npy').astype(np.int64)
        r0 = np.argmax(n0[n0] == np.arange(len(n0))[:, None, None], axis=2)
        d0 = p0[n0] - p0[:, None, :]; d0 -= L0 * np.round(d0 / L0)
    else:
        p0, n0, r0, d0, L0 = load_net({'glass': 'glass', 'diamond': 'diamond', 'd216': 'diamond_small'}[name])
    bip0 = bipartite_sign(n0) if name in ('diamond', 'd216') else None
    rev0 = np.ascontiguousarray(r0, dtype=np.int32)
    if T == 1:
        P, NB, LB = p0 % L0, n0.astype(np.int32), L0
    else:
        P, NB, LB = lean.tile(p0, n0, rev0, L0, T)
    bip = None if bip0 is None else np.tile(bip0, T ** 3)
    del n0, r0; gc.collect()
    return P, NB, rev0, np.ascontiguousarray(d0), LB, L0, bip

def gridn(LB): return int(max(32, 16 * round(LB / 5.0866 / 16)))

def one(pos, nb, rev0, d0, LB, bip, c, vf, dn, SIG, c0, Tn, log):
    u = D2[dn]; g = 1 / np.sqrt(1 - vf ** 2); kv = g * eta * vf / c * u; kn = np.linalg.norm(kv)
    w = omega(d0, kv, eta); dk = 1e-4 * u; vg = (omega(d0, kv + dk, eta) - omega(d0, kv - dk, eta)) / 2e-4
    t0 = time.time()
    phi, e_in, e_out, info = lean.build_packet(pos, nb, d0, LB, c0, kv, w, SIG, eta, mu_of_k, bip)
    tb = time.time() - t0; log(f"built: cg {info} e_in {e_in:.1e} e_out {e_out:.1e} {tb:.0f}s rss {rss():.2f} GB peak {peak():.2f} GB")
    n = gridn(LB)
    shell, fwd = masks(n, LB, kv); grid = np.zeros((n, n, n), np.complex128)
    Wt, wmin, wmax = H2.tracker_weight(fwd, n, LB, d0)
    out = np.zeros_like(phi); ce = np.cos(eta) / 2; em = np.exp(-1j * eta); z3 = np.zeros(3); rc2 = (0.45 * LB) ** 2
    cen = c0.copy(); un = c0.copy(); cc = c0.copy(); uc = c0.copy(); co = c0.copy(); uo = c0.copy()
    rec = {k: [] for k in ['ts', 'xs', 'xc', 'xo', 'ph', 'nm', 'rms', 'rms_c', 's', 'off', 'fwd_w']}
    tstep = 0.0; tmeas = 0.0
    for t in range(Tn + 1):
        if t % 5 == 0:
            tm = time.time()
            for _ in range(3):
                m = moments(phi, pos, cen, LB, rc2, kv, z3); sh = m[1:4] / m[0]; cen = (cen + sh) % LB; un = un + sh
            m = moments(phi, pos, cen, LB, rc2, kv, un - c0)
            grid[:] = 0; cic(phi, pos, LB, n, grid); F = np.fft.fftn(grid); P = np.abs(F) ** 2
            Ptot = P.sum(); Psh = P[shell].sum(); Pf = P[fwd].sum()
            rho_c = np.abs(np.fft.ifftn(F * Wt)) ** 2; rho_o = np.abs(np.fft.ifftn(np.where(fwd, F, 0))) ** 2
            for _ in range(3):
                gm = grid_centroid(rho_c, cc, LB); sh = gm[1:4] / gm[0]; cc = (cc + sh) % LB; uc = uc + sh
                go = grid_centroid(rho_o, co, LB); sh = go[1:4] / go[0]; co = (co + sh) % LB; uo = uo + sh
            gm = grid_centroid(rho_c, cc, LB)
            rec['ts'].append(t); rec['xs'].append(float((un - c0) @ (-u))); rec['xc'].append(float((uc - c0) @ (-u))); rec['xo'].append(float((uo - c0) @ (-u)))
            rec['ph'].append(float(np.angle(m[9] + 1j * m[10]))); rec['nm'].append(float(m[0])); rec['rms'].append(float(np.sqrt(m[8] / m[0])))
            rec['rms_c'].append(float(np.sqrt(gm[4] / gm[0]))); rec['s'].append(float(1 - Pf / Psh)); rec['off'].append(float(1 - Psh / Ptot)); rec['fwd_w'].append(float(Pf / Ptot))
            tmeas += time.time() - tm
            if t % 100 == 0: log(f"tick {t}: ds {100*(rec['s'][-1]-rec['s'][0]):+.3f} pp, rss {rss():.2f} GB")
        if t < Tn:
            ts_ = time.time(); lean.step(phi, nb, rev0, ce, em, out); phi, out = out, phi; tstep += time.time() - ts_
    del phi, out; gc.collect()
    rec.update(kn=float(kn), e_in=e_in, e_out=e_out, cg_info=int(info), t_build=tb, t_total=time.time() - t0, grid_n=n, w_min=wmin, w_max=wmax,
               vg_carrier=float(vg / c), t_step=tstep / max(Tn, 1), t_meas=tmeas)
    return rec

if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--net', required=True); ap.add_argument('--tile', type=int, default=2)
    ap.add_argument('--point', type=int, default=1); ap.add_argument('--ticks', type=int, default=600); ap.add_argument('--dir', default='100')
    ap.add_argument('--start', type=int, default=0); ap.add_argument('--out'); a = ap.parse_args()
    log = lambda s: print(f"[{time.strftime('%H:%M:%S')}] {s}", flush=True)
    t0 = time.time(); P, NB, rev0, d0, LB, L0, bip = load(a.net, a.tile)
    log(f"{a.net} x{a.tile}: N={len(NB):,} box={LB:.1f} A loaded {time.time()-t0:.0f}s rss {rss():.2f} GB peak {peak():.2f} GB")
    c = CT.get({'glass': 'glass', 'diamond': 'diamond', 'g216': 'raw_tiled_test', 'd216': 'diamond_small'}[a.net], 1.3562)
    pt = POINTS[a.point - 1]
    c0 = np.array(STARTS[a.start]) * L0          # starts in units of the base box, as in 4f-4i
    Tn = a.ticks; ts = min(TS, Tn // 2)
    r = one(P, NB, rev0, d0, LB, bip, c, pt['vf'], a.dir, pt['sig'], c0, Tn, log)
    sc = H2.score(r, c, Tn, ts); sc.update(net=a.net, tile=a.tile, point=a.point, dir=a.dir, N=len(NB), box=LB, vf=pt['vf'], sig=pt['sig'])
    sc['inv_tau'] = decay(r, min(100, Tn // 4)); vp = sc['coh']['v'] * c; sc['inv_ell'] = float(sc['inv_tau'] / vp)
    sc['peak_GB'] = peak(); sc['t_step'] = r['t_step']; sc['t_meas'] = r['t_meas']
    log(f"RESULT particle-only v/c {sc['coh']['v']:.5f} whole-box {sc['whole']['v']:.5f} | 1/tau {sc['inv_tau']:.3e} Gamma {sc['Gamma']:+.3f} max|ds| {sc['ds_max']:.3f} s0 {sc['s0']:.4f} | "
        f"clock {sc['clock_err']:+.2f}% norm {sc['norm']:.1e} | build {sc['t_build']:.0f}s step {sc['t_step']:.3f}s/tick total {sc['t_total']:.0f}s | PEAK {sc['peak_GB']:.2f} GB")
    if a.out: os.makedirs(os.path.dirname(a.out) or '.', exist_ok=True); json.dump(sc, open(a.out, 'w'))
