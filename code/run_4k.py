"""Experiment 4k: direction spread on untiled glass s5 vs tiled s3/s4 (6 directions, points 1 and 5, 2 starts).
Usage: python run_4k.py JOB_INDEX [--out DIR]   |   python run_4k.py --list"""
import numpy as np, json, time, os, argparse
import run_4g
run_4g.D2.update({'010': np.array([0., 1, 0]), '001': np.array([0., 0, 1]),
                  '312': np.array([3., 1, 2]) / np.sqrt(14), '231': np.array([2., 3, 1]) / np.sqrt(14)})
import run_4h2 as H2
import run_4i as I
from run_4i_rep import points_for
from run_4g import STARTS
from run_4h import CT, T, TS

DIRS6 = ['100', '010', '001', '123', '312', '231']
def jobs():
    J = []
    P5 = [p for p in points_for(CT['glass5']) if p['point'] in (1, 5)]
    for p in P5:
        for dn in DIRS6:
            for s in STARTS[:2]: J.append(dict(net='glass5', dir=dn, start=s, **p))
    for net in ['glass', 'glass4']:
        for p in [q for q in points_for(CT[net]) if q['point'] in (1, 5)]:
            for dn in ['010', '001', '312', '231']:
                for s in STARTS[:2]: J.append(dict(net=net, dir=dn, start=s, **p))
    return J

if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('job', nargs='?', type=int); ap.add_argument('--out', default='../results/exp4k'); ap.add_argument('--list', action='store_true'); a = ap.parse_args()
    if a.list:
        for i, j in enumerate(jobs()): print(i, j)
        raise SystemExit
    spec = jobs()[a.job]
    t0 = time.time(); pos, nb, rev, d0, LB, LB0, bip = H2.load(spec['net']); A = H2.adj(nb); c = CT[spec['net']]
    print(f"job {a.job} {spec}: c={c} N={len(nb)} box={LB:.1f} base={LB0:.1f} loaded {time.time()-t0:.0f}s", flush=True)
    c0 = np.array(spec['start']) * LB0
    r = H2.one(pos, nb, rev, d0, LB, A, bip, c, spec['vf'], spec['dir'], spec['sig'], c0, T)
    sc = H2.score(r, c, T, TS); sc.update(spec); sc['N'] = len(nb); sc['box'] = LB; sc['c'] = c
    sc['inv_tau'] = I.decay(r, 100); vp = sc['coh']['v'] * c; sc['inv_ell'] = float(sc['inv_tau'] / vp); sc['k0'] = float(r['kn'])
    print(f"{spec['net']} point {spec['point']} dir {spec['dir']} launch {spec['vf']:.5f}: particle-only v/c {sc['coh']['v']:.5f} (whole-box {sc['whole']['v']:.5f}) "
          f"1/tau {sc['inv_tau']:.4e} 1/ell {sc['inv_ell']:.4e} | max|ds| {sc['ds_max']:.2f} s0 {sc['s0']:.4f} | clock {sc['clock_err']:+.2f}% norm {sc['norm']:.1e} | total {sc['t_total']:.0f}s", flush=True)
    os.makedirs(a.out, exist_ok=True); json.dump(sc, open(f"{a.out}/job_{a.job:02d}.json", 'w'))
