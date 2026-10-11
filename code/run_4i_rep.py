"""Experiment 4i-R: replication of 4i on the independent glass s4 (same k0 and sigma as 4i point by point).
Usage: python run_4i_rep.py JOB_INDEX [--out DIR]"""
import numpy as np, json, time, os, argparse, sys
import run_4h2 as H2
import run_4i as I
from run_4g import STARTS, eta
from run_4h import CT, T, TS

C3 = CT['glass']

def points_for(c):
    """Same k0 and sigma as 4i (defined with c_s3); launch speed re-derived for speed limit c."""
    P = []
    for p in I.POINTS:
        k0 = (p['vf'] / np.sqrt(1 - p['vf'] ** 2)) * eta / C3
        gv = k0 * c / eta
        P.append(dict(point=p['point'], vf=float(gv / np.sqrt(1 + gv ** 2)), sig=p['sig'], k0_target=float(k0)))
    return P

def jobs():
    c4 = CT['glass4']; J = []
    for p in points_for(c4):
        for dn in ['100', '123']:
            for s in STARTS[:3]: J.append(dict(net='glass4', dir=dn, start=s, set='G', **p))
    p5 = I.POINTS[4]; J.append(dict(net='diamond', dir='100', start=STARTS[0], set='C', **p5))   # pipeline check = 4i job 34
    return J

if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('job', type=int); ap.add_argument('--out', default='../results/exp4i_rep'); a = ap.parse_args()
    spec = jobs()[a.job]
    t0 = time.time(); pos, nb, rev, d0, LB, LB0, bip = H2.load(spec['net']); A = H2.adj(nb); c = CT[spec['net']]
    print(f"job {a.job} {spec}: c={c} N={len(nb)} box={LB:.1f} loaded {time.time()-t0:.0f}s", flush=True)
    c0 = np.array(spec['start']) * LB0
    r = H2.one(pos, nb, rev, d0, LB, A, bip, c, spec['vf'], spec['dir'], spec['sig'], c0, T)
    sc = H2.score(r, c, T, TS); sc.update(spec); sc['N'] = len(nb); sc['box'] = LB; sc['c'] = c
    sc['inv_tau'] = I.decay(r, 100); vp = sc['coh']['v'] * c
    sc['ell'] = float(vp / sc['inv_tau']) if sc['inv_tau'] > 0 else None; sc['inv_ell'] = float(sc['inv_tau'] / vp); sc['k0'] = float(r['kn'])
    print(f"point {spec['point']} launch {spec['vf']:.5f} sigma {spec['sig']:.1f} k0 {sc['k0']:.6f}: particle-only v/c {sc['coh']['v']:.5f} (whole-box {sc['whole']['v']:.5f}) "
          f"1/tau {sc['inv_tau']:.4e} | Gamma {sc['Gamma']:+.3f} max|ds| {sc['ds_max']:.2f} s0 {sc['s0']:.4f} | clock {sc['clock_err']:+.2f}% norm {sc['norm']:.1e} | total {sc['t_total']:.0f}s", flush=True)
    os.makedirs(a.out, exist_ok=True); json.dump(sc, open(f"{a.out}/job_{a.job:02d}.json", 'w'))
