"""Experiment 4i: glass scattering vs wavelength with scale-matched packets (sigma ∝ wavelength).
Usage: python run_4i.py JOB_INDEX [--out DIR]   or   --test NET --point P --T T"""
import numpy as np, json, time, os, argparse
import run_4h2 as H2
from run_4g import STARTS, eta
from run_4h import CT, T, TS

R4 = 2 ** 0.25
POINTS = []
for j in range(5):
    gv = (4 / 3) / R4 ** j
    POINTS.append(dict(point=j + 1, vf=float(gv / np.sqrt(1 + gv ** 2)), sig=float(40.0 * R4 ** j)))

def jobs():
    J = []
    for p in POINTS:
        for dn in ['100', '123']:
            for s in STARTS[:3]: J.append(dict(net='glass', dir=dn, start=s, set='G', **p))
    for p in POINTS: J.append(dict(net='diamond', dir='100', start=STARTS[0], set='C', **p))
    return J

def decay(r, t0=100):
    t = np.array(r['ts'], float); s = np.array(r['s']); f = 1 - (s - s[0]); m = t >= t0
    return float(np.polyfit(t[m], -np.log(f[m]), 1)[0])

if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('job', nargs='?', type=int); ap.add_argument('--out', default='../results/exp4i')
    ap.add_argument('--test'); ap.add_argument('--point', type=int, default=1); ap.add_argument('--T', type=int, default=T); ap.add_argument('--dir', default='100')
    a = ap.parse_args()
    spec = dict(net=a.test, dir=a.dir, start=STARTS[1], set='test', **POINTS[a.point - 1]) if a.test else jobs()[a.job]
    t0 = time.time(); pos, nb, rev, d0, LB, LB0, bip = H2.load(spec['net']); A = H2.adj(nb); c = CT.get(spec['net'], 1.3562)
    print(f"job {a.job} {spec}: N={len(nb)} box={LB:.1f} loaded {time.time()-t0:.0f}s", flush=True)
    c0 = np.array(spec['start']) * LB0; Tn = a.T if a.test else T; ts = min(TS, Tn // 2)
    r = H2.one(pos, nb, rev, d0, LB, A, bip, c, spec['vf'], spec['dir'], spec['sig'], c0, Tn)
    sc = H2.score(r, c, Tn, ts); sc.update(spec); sc['N'] = len(nb); sc['box'] = LB
    sc['inv_tau'] = decay(r, min(100, Tn // 4)); vp = sc['coh']['v'] * c
    sc['ell'] = float(vp / sc['inv_tau']) if sc['inv_tau'] > 0 else None; sc['inv_ell'] = float(sc['inv_tau'] / vp); sc['k0'] = float(r['kn'])
    print(f"point {spec['point']} v/c launch {spec['vf']:.4f} sigma {spec['sig']:.1f}: particle-only v/c {sc['coh']['v']:.4f} (whole-box {sc['whole']['v']:.4f}) "
          f"1/tau {sc['inv_tau']:.3e} ell {sc['ell'] if sc['ell'] is None else round(sc['ell'])} A | Gamma {sc['Gamma']:+.3f} max|ds| {sc['ds_max']:.2f} s0 {sc['s0']:.4f} | "
          f"clock {sc['clock_err']:+.2f}% norm {sc['norm']:.1e} build {sc['t_build']:.0f}s total {sc['t_total']:.0f}s", flush=True)
    if not a.test:
        os.makedirs(a.out, exist_ok=True); json.dump(sc, open(f"{a.out}/job_{a.job:02d}.json", 'w'))
    else:
        json.dump(sc, open(f"test_4i_{a.test}_p{a.point}_{a.dir}.json", 'w'))
