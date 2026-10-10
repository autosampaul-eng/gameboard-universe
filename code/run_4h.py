"""Experiment 4h: separate the coherent particle from scattered waves using the wave's direction spectrum.
Usage: python run_4h.py JOB_INDEX [--out DIR]     (or --test NAME for method checks on 216k networks)"""
import numpy as np, json, time, sys, os, gc, argparse
from numba import njit, prange
from exp4 import load_net, omega, step
from run_4g import load_tiled, adj, build_packet, moments, D2, STARTS, eta

T = 600; TS = 250; TH = 425; CONE = np.cos(np.radians(30.0)); SHELL = 0.5
CT = {'glass': 1.5495, 'diamond': 1.3562, 'glass_small_box': 1.5495, 'raw_tiled_test': 1.529, 'diamond_small': 1.3562}

def jobs():
    J = []
    for dn in ['100', '123']:
        for s in STARTS: J.append(dict(net='glass', vf=0.8, dir=dn, start=s, set='A'))
    for dn in ['100', '123']:
        for s in STARTS[:3]: J.append(dict(net='glass', vf=0.5, dir=dn, start=s, set='B'))
    for vf in [0.8, 0.5]:
        for dn in ['100', '123']: J.append(dict(net='diamond', vf=vf, dir=dn, start=STARTS[0], set='C'))
    for s in STARTS[:3]: J.append(dict(net='glass_small_box', vf=0.8, dir='100', start=s, set='D'))
    return J

def load(name):
    if name in ('glass', 'diamond', 'raw_tiled_test'):
        return load_tiled(name)
    base = {'glass_small_box': 'glass', 'diamond_small': 'diamond_small'}[name]
    pos, nb, rev, d, LB = load_net(base)
    from exp4 import bipartite_sign
    bip = bipartite_sign(nb) if 'diamond' in name else None
    return pos, nb.astype(np.int32), rev.astype(np.int32), d, LB, LB, bip

@njit(cache=True)
def cic(phi, pos, LB, n, out):
    h = LB / n
    for i in range(phi.shape[0]):
        S = phi[i, 0] + phi[i, 1] + phi[i, 2] + phi[i, 3]
        x = pos[i, 0] / h; y = pos[i, 1] / h; z = pos[i, 2] / h
        i0 = int(np.floor(x)); j0 = int(np.floor(y)); k0 = int(np.floor(z))
        fx = x - i0; fy = y - j0; fz = z - k0
        for a in range(2):
            wx = fx if a else 1 - fx; ia = (i0 + a) % n
            for b in range(2):
                wy = fy if b else 1 - fy; jb = (j0 + b) % n
                for c in range(2):
                    wz = fz if c else 1 - fz; kc = (k0 + c) % n
                    out[ia, jb, kc] += S * (wx * wy * wz)

@njit(parallel=True, cache=True)
def grid_centroid(rho, cen, LB):
    n = rho.shape[0]; h = LB / n
    w = 0.0; sx = 0.0; sy = 0.0; sz = 0.0; q = 0.0
    for i in prange(n):
        for j in range(n):
            for k in range(n):
                r = rho[i, j, k]
                dx = i * h - cen[0]; dx -= LB * np.round(dx / LB)
                dy = j * h - cen[1]; dy -= LB * np.round(dy / LB)
                dz = k * h - cen[2]; dz -= LB * np.round(dz / LB)
                w += r; sx += r * dx; sy += r * dy; sz += r * dz; q += r * (dx * dx + dy * dy + dz * dz)
    return np.array([w, sx, sy, sz, q])

def masks(n, LB, kv):
    k1 = 2 * np.pi * np.fft.fftfreq(n, d=LB / n)
    KX, KY, KZ = np.meshgrid(k1, k1, k1, indexing='ij'); K = np.sqrt(KX ** 2 + KY ** 2 + KZ ** 2)
    k0 = np.linalg.norm(kv); u = kv / k0
    shell = np.abs(K - k0) < SHELL * k0
    cos = (KX * u[0] + KY * u[1] + KZ * u[2]) / np.maximum(K, 1e-30)
    fwd = shell & (cos > CONE)
    return shell, fwd

def one(pos, nb, rev, d0, LB, A, bip, c, vf, dn, SIG, c0, Tn=T):
    u = D2[dn]; g = 1 / np.sqrt(1 - vf ** 2); kv = g * eta * vf / c * u; kn = np.linalg.norm(kv)
    w = omega(d0, kv, eta); dk = 1e-4 * u; vg = (omega(d0, kv + dk, eta) - omega(d0, kv - dk, eta)) / 2e-4
    t0 = time.time()
    phi, e_in, e_out, info = build_packet(pos, nb, d0, LB, A, c0, kv, -vg * u, w, SIG, bip)
    tb = time.time() - t0
    n = 128 if LB > 500 else (64 if LB > 250 else 32)
    shell, fwd = masks(n, LB, kv); grid = np.zeros((n, n, n), np.complex128)
    out = np.zeros_like(phi); ce = np.cos(eta) / 2; em = np.exp(-1j * eta); z3 = np.zeros(3); rc2 = (0.45 * LB) ** 2
    cen = c0.copy(); un = c0.copy(); cc = c0.copy(); uc = c0.copy()
    rec = {k: [] for k in ['ts', 'xs', 'xc', 'ph', 'nm', 'rms', 'rms_c', 's', 'off', 'fwd_w']}
    for t in range(Tn + 1):
        if t % 5 == 0:
            for _ in range(3):
                m = moments(phi, pos, cen, LB, rc2, kv, z3); sh = m[1:4] / m[0]; cen = (cen + sh) % LB; un = un + sh
            m = moments(phi, pos, cen, LB, rc2, kv, un - c0)
            grid[:] = 0; cic(phi, pos, LB, n, grid); F = np.fft.fftn(grid); P = np.abs(F) ** 2
            Ptot = P.sum(); Psh = P[shell].sum(); Pf = P[fwd].sum()
            rho_c = np.abs(np.fft.ifftn(np.where(fwd, F, 0))) ** 2
            for _ in range(3):
                gm = grid_centroid(rho_c, cc, LB); sh = gm[1:4] / gm[0]; cc = (cc + sh) % LB; uc = uc + sh
            gm = grid_centroid(rho_c, cc, LB)
            rec['ts'].append(t); rec['xs'].append(float((un - c0) @ (-u))); rec['xc'].append(float((uc - c0) @ (-u)))
            rec['ph'].append(float(np.angle(m[9] + 1j * m[10]))); rec['nm'].append(float(m[0])); rec['rms'].append(float(np.sqrt(m[8] / m[0])))
            rec['rms_c'].append(float(np.sqrt(gm[4] / gm[0]))); rec['s'].append(float(1 - Pf / Psh)); rec['off'].append(float(1 - Psh / Ptot)); rec['fwd_w'].append(float(Pf / Ptot))
        if t < Tn: step(phi, nb, rev, ce, em, out); phi, out = out, phi
    del phi, out; gc.collect()
    rec.update(kn=float(kn), e_in=e_in, e_out=e_out, cg_info=int(info), t_build=tb, t_total=time.time() - t0, grid_n=n)
    return rec

def speed_stats(t, x, c, Tn=T, ts=TS, th=TH):
    m = t >= ts; v = np.polyfit(t[m], x[m], 1)[0]
    h1 = (t >= ts) & (t <= th); h2 = t >= th
    v1 = np.polyfit(t[h1], x[h1], 1)[0]; v2 = np.polyfit(t[h2], x[h2], 1)[0]
    wc = []; wv = []
    for a0 in range(0, Tn - 39, 20):
        mm = (t >= a0) & (t <= a0 + 40); wv.append(np.polyfit(t[mm], x[mm], 1)[0]); wc.append(a0 + 20)
    wc = np.array(wc, float); wv = np.array(wv); ms = wc - 20 >= ts
    p = np.polyfit(wc[ms], wv[ms], 1)
    return dict(v=float(v / c), half=float(100 * (v2 / v1 - 1)), D=float(100 * p[0] * 100 / v),
                W=float(100 * np.sqrt(((wv[ms] - np.polyval(p, wc[ms])) ** 2).mean()) / v), win_t=wc.tolist(), win_v=(wv / c).tolist()), v

def score(r, c, Tn=T, ts=TS):
    t = np.array(r['ts'], float)
    whole, _ = speed_stats(t, np.array(r['xs']), c, Tn, ts, (ts + Tn) / 2)
    coh, vc = speed_stats(t, np.array(r['xc']), c, Tn, ts, (ts + Tn) / 2)
    m = t >= ts
    wm = np.polyfit(t[m], np.unwrap(r['ph'])[m], 1)[0]; clock = wm - r['kn'] * vc
    s = np.array(r['s']); ds = s - s[0]
    G = 100 * np.polyfit(t[m], 100 * ds[m], 1)[0]           # percentage points per 100 ticks
    ds_ts = float(100 * np.interp(ts, t, ds))
    mfp = float(vc * 100 / (G / 100)) if G > 0 else None
    return dict(whole=whole, coh=coh, clock_err=float(100 * (clock / eta / np.sqrt(1 - (vc / c) ** 2) - 1)), Gamma=float(G), ds_startup=ds_ts,
                ds_max=float(100 * np.abs(ds).max()), s0=float(s[0]), mfp=mfp, norm=float(abs(r['nm'][-1] / r['nm'][0] - 1)),
                rms_end=r['rms'][-1], rms_c_end=r['rms_c'][-1], off_end=r['off'][-1], off0=r['off'][0], ds_series=(100 * ds).tolist(), ts=r['ts'],
                e_out=r['e_out'], cg_info=r['cg_info'], t_build=r['t_build'], t_total=r['t_total'], grid_n=r['grid_n'])

if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('job', nargs='?', type=int); ap.add_argument('--out', default='../results/exp4h')
    ap.add_argument('--test'); ap.add_argument('--T', type=int, default=T); ap.add_argument('--sig', type=float, default=40.0); ap.add_argument('--vf', type=float, default=0.8)
    ap.add_argument('--dir', default='100'); a = ap.parse_args()
    if a.test:
        spec = dict(net=a.test, vf=a.vf, dir=a.dir, start=STARTS[1], set='test')
    else:
        spec = jobs()[a.job]
    t0 = time.time(); pos, nb, rev, d0, LB, LB0, bip = load(spec['net']); A = adj(nb); c = CT[spec['net']]
    print(f"job {a.job} {spec}: N={len(nb)} box={LB:.1f} loaded {time.time()-t0:.0f}s", flush=True)
    c0 = np.array(spec['start']) * LB0
    Tn = a.T if a.test else T; ts = min(TS, Tn // 2)
    r = one(pos, nb, rev, d0, LB, A, bip, c, spec['vf'], spec['dir'], a.sig if a.test else 40.0, c0, Tn)
    sc = score(r, c, Tn, ts); sc.update(spec); sc['N'] = len(nb); sc['box'] = LB
    print(f"particle-only v/c {sc['coh']['v']:.4f} (whole-box {sc['whole']['v']:.4f}) D {sc['coh']['D']:+.2f} halves {sc['coh']['half']:+.2f}% | "
          f"s0 {sc['s0']:.4f} startup ds {sc['ds_startup']:+.2f}pp Gamma {sc['Gamma']:+.3f}pp/100t max|ds| {sc['ds_max']:.2f} | clock {sc['clock_err']:+.2f}% "
          f"norm {sc['norm']:.1e} rms {sc['rms_end']:.0f}/{sc['rms_c_end']:.0f} build {sc['t_build']:.0f}s total {sc['t_total']:.0f}s", flush=True)
    if not a.test:
        os.makedirs(a.out, exist_ok=True); json.dump(sc, open(f"{a.out}/job_{a.job:02d}.json", 'w'))
    else:
        json.dump(sc, open(f"test_4h_{a.test}_{a.dir}_{a.vf}.json", 'w'))
