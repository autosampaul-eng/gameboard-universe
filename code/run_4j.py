"""Experiment 4j: separating wavelength from speed (eta varied). 3 speeds x 7 wavelengths, 3x3x3 tiled box.
Usage: python run_4j.py JOB [--out DIR] | --list | --table | --test NET JOB_TEMPLATE --tile T --ticks N"""
import numpy as np, json, time, os, argparse
import run_4g
run_4g.D2.update({'312': np.array([3., 1, 2]) / np.sqrt(14)})
import run_big as RB, run_4h2 as H2
from run_4i import POINTS, decay
from run_4g import STARTS
from run_4h import CT, TS

R = 2 ** 0.25; ETA0 = 0.15; T = 600; TILE = 3
SPEEDS = [0, 1, 2]                      # speed index i: gamma*beta = (4/3) / R^(2i)  ->  4i points 1, 3, 5
WAVES = list(range(7))                  # wavelength index j: k0 = (4/3) * 0.15 / (c_s3 * R^j)
DIRS = [('100', 0), ('123', 1), ('312', 2)]
C = {'glass': CT['glass'], 'glass4': CT['glass4'], 'diamond': CT['diamond']}

def cell(i, j, net):
    vf = POINTS[2 * i]['vf']; gb = vf / np.sqrt(1 - vf ** 2)
    eta = ETA0 * R ** (2 * i - j)                       # s3 (and diamond): exactly 0.15 when j == 2i
    if net == 'glass4': eta *= C['glass4'] / C['glass']  # same k0 as s3 at the same speed
    sig = 40.0 * R ** j
    return dict(i=i, j=j, vf=vf, gb=gb, eta=eta, sig=sig, k0=gb * eta / C[net], lam_bonds=2 * np.pi / (gb * eta / C[net]) / 2.735)

def v_pkt_analytic(gb, eta, c, sig, nq=48):
    """Registered reference speed (units of c): continuum group velocity c^2 kx / sqrt(eta^2 + c^2 k^2) averaged over the
    packet's wavenumber distribution, a 3D Gaussian of std 1/(2 sigma) per component about k0 = gb*eta/c along x. Gauss-Hermite."""
    x, w = np.polynomial.hermite_e.hermegauss(nq); w = w / w.sum(); sk = 1 / (2 * sig); k0 = gb * eta / c
    KX, KY, KZ = np.meshgrid(k0 + sk * x, sk * x, sk * x, indexing='ij'); W = w[:, None, None] * w[None, :, None] * w[None, None, :]
    return float((W * c * KX / np.sqrt(eta ** 2 + c * c * (KX ** 2 + KY ** 2 + KZ ** 2))).sum())

def jobs():
    J = []
    for i in SPEEDS:
        for j in WAVES:
            for net in ('glass', 'glass4'):
                for dn, st in DIRS: J.append(dict(net=net, dir=dn, start=st, set='G', **cell(i, j, net)))
    for i in SPEEDS:
        for j in WAVES: J.append(dict(net='diamond', dir='100', start=0, set='C', **cell(i, j, 'diamond')))
    J.append(dict(net='diamond', dir='100', start=0, set='REG', **cell(2, 4, 'diamond')))   # unchanged default code path
    for q in J: q['v_pkt_an'] = v_pkt_analytic(q['gb'], q['eta'], C[q['net']], q['sig'])
    return J

def run(spec, tile=TILE, Tn=T, log=print):
    t0 = time.time(); P, NB, rev0, d0, LB, L0, bip = RB.load(spec['net'], tile); c = C.get(spec['net'], spec.get('c'))
    log(f"{spec}: N={len(NB):,} box={LB:.1f} loaded {time.time()-t0:.0f}s")
    c0 = np.array(STARTS[spec['start']]) * L0
    kw = {} if spec['set'] == 'REG' else dict(eta_=spec['eta'])
    r = RB.one(P, NB, rev0, d0, LB, bip, c, spec['vf'], spec['dir'], spec['sig'], c0, Tn, log, **kw)
    ts = min(TS, Tn // 2); sc = H2.score(r, c, Tn, ts); sc.update(spec); sc['N'] = len(NB); sc['box'] = LB; sc['c'] = c
    sc['clock_err'] = float(100 * ((sc['clock_err'] / 100 + 1) * ETA0 / spec['eta'] - 1))   # score() divides by the default eta
    sc['inv_tau'] = decay(r, min(100, Tn // 4)); vp = sc['coh']['v'] * c; sc['inv_ell'] = float(sc['inv_tau'] / vp)
    sc['k0_run'] = float(r['kn']); sc['v_pkt_spec'] = r['v_pkt']
    sc['delta'] = float(1 - sc['coh']['v'] / spec['v_pkt_an']); sc['delta_spec'] = float(1 - sc['coh']['v'] / r['v_pkt']); sc['shortfall_raw'] = float(1 - sc['coh']['v'] / spec['vf'])
    sc['t_step'] = r['t_step']; sc['peak_GB'] = RB.peak()
    log(f"RESULT {spec['net']} i{spec['i']} j{spec['j']} {spec['dir']} eta {spec['eta']:.4f} k0 {sc['k0_run']:.5f}: v_p {sc['coh']['v']:.5f} (launch {spec['vf']:.5f}, "
        f"packet {spec['v_pkt_an']:.5f}) delta {100*sc['delta']:+.3f}% raw {100*sc['shortfall_raw']:+.3f}% | 1/ell {sc['inv_ell']:.4e} | max|ds| {sc['ds_max']:.3f} "
        f"1/tau {sc['inv_tau']:.2e} whole {sc['whole']['v']:.5f} | clock {sc['clock_err']:+.2f}% norm {sc['norm']:.1e} | {sc['t_total']:.0f}s peak {sc['peak_GB']:.2f} GB")
    return sc

if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('job', nargs='?', type=int); ap.add_argument('--out', default='../results/exp4j')
    ap.add_argument('--list', action='store_true'); ap.add_argument('--table', action='store_true')
    ap.add_argument('--test'); ap.add_argument('--tile', type=int, default=TILE); ap.add_argument('--ticks', type=int, default=T); a = ap.parse_args()
    J = jobs()
    if a.list:
        for n, q in enumerate(J): print(n, q['set'], q['net'], q['i'], q['j'], q['dir'], q['start'], f"eta {q['eta']:.5f} vf {q['vf']:.5f} k0 {q['k0']:.5f} sig {q['sig']:.1f} vpkt {q['v_pkt_an']:.5f}")
        raise SystemExit
    if a.table:
        print('| Speed | j | Wavelength (bonds) | Packet speed / launch, s3 | s4 | diamond |'); print('|---|---|---|---|---|---|')
        for i in SPEEDS:
            for j in WAVES:
                q = {net: [x for x in J if x['net'] == net and x['i'] == i and x['j'] == j and x['set'] != 'REG'][0] for net in ('glass', 'glass4', 'diamond')}
                print(f"| {q['glass']['vf']:.3f}c | {j+1} | {q['glass']['lam_bonds']:.1f} | " + ' | '.join(f"{100*(q[n]['v_pkt_an']/q[n]['vf']-1):+.3f}%" for n in ('glass', 'glass4', 'diamond')) + ' |')
        raise SystemExit
    log = lambda s: print(f"[{time.strftime('%H:%M:%S')}] {s}", flush=True)
    if a.test:
        spec = dict(J[a.job]); spec['net'] = a.test; spec['c'] = {'g216': 1.529, 'd216': 1.3562}[a.test]; C[a.test] = spec['c']
        sc = run(spec, a.tile, a.ticks, log); raise SystemExit
    sc = run(J[a.job], log=log)
    os.makedirs(a.out, exist_ok=True); json.dump(sc, open(f"{a.out}/job_{a.job:03d}.json", 'w'))
