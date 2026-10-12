"""Exp 4k scoring, exactly as pre-registered in DESIGN_exp4k.md."""
import json, glob, numpy as np, sys
from run_4g import STARTS
S2 = [list(s) for s in STARTS[:2]]
def L(p): return [json.load(open(f)) for f in sorted(glob.glob(p))]
K = L('res_4k/job_*.json'); I3 = L('res_4i/job_*.json'); I4 = L('res_4i_rep/job_*.json')
DIRS = ['100', '010', '001', '123', '312', '231']
def runs(net, point):
    out = {d: [] for d in DIRS}
    src = [r for r in K if r['net'] == net and r['point'] == point]
    if net != 'glass5':
        old = I3 if net == 'glass' else I4
        src += [r for r in old if r.get('set') == 'G' and r['point'] == point and r['dir'] in ('100', '123') and list(r['start']) in S2]
    for r in src: out[r['dir']].append(r)
    return out
def spread(R):
    m = np.array([np.mean([r['inv_ell'] for r in R[d]]) for d in DIRS])
    se2 = np.array([np.var([r['inv_ell'] for r in R[d]], ddof=1) / len(R[d]) for d in DIRS])
    return float(np.sqrt(max(m.var(ddof=1) - se2.mean(), 0)) / m.mean()), m, float(m.mean())
def boot(R, nb=4000, rng=np.random.default_rng(0)):
    out = []
    for _ in range(nb):
        RR = {d: [R[d][i] for i in rng.integers(0, len(R[d]), len(R[d]))] for d in DIRS}
        out.append(spread(RR)[0])
    return np.array(out)
res = {}
for point in (5, 1):
    row = {}
    for net in ('glass5', 'glass', 'glass4'):
        R = runs(net, point); n = {d: len(R[d]) for d in DIRS}; assert all(v == 2 for v in n.values()), (net, point, n)
        A, m, mean = spread(R); row[net] = dict(A=A, m=m.tolist(), mean=mean, boot=boot(R))
        clocks = [r['clock_err'] for d in DIRS for r in R[d]]; norms = [r['norm'] for d in DIRS for r in R[d]]
        row[net]['clock_range'] = [min(clocks), max(clocks)]; row[net]['norm_max'] = max(norms)
        row[net]['v'] = float(np.mean([r['coh']['v'] for d in DIRS for r in R[d]])); row[net]['vf'] = R['100'][0]['vf']
    At = np.sqrt((row['glass']['A'] ** 2 + row['glass4']['A'] ** 2) / 2); Rr = row['glass5']['A'] / At if At > 0 else np.inf
    bR = row['glass5']['boot'] / np.sqrt((row['glass']['boot'] ** 2 + row['glass4']['boot'] ** 2) / 2)
    verdict = 'shrinks (finite-sample)' if Rr <= 0.6 else ('does not shrink' if Rr >= 0.85 else 'inconclusive')
    tiled_mean = (row['glass']['mean'] + row['glass4']['mean']) / 2; Sratio = row['glass5']['mean'] / tiled_mean
    print(f"\n=== point {point} ===")
    for net in ('glass5', 'glass', 'glass4'):
        r = row[net]; print(f"{net:7s} A {100*r['A']:.1f}% (boot 68%: {100*np.percentile(r['boot'],16):.1f}-{100*np.percentile(r['boot'],84):.1f}) mean 1/ell {r['mean']:.4e} | per dir (vs mean %): "
              + ' '.join(f"{d}:{100*(x/r['mean']-1):+.1f}" for d, x in zip(DIRS, r['m'])) + f" | v_p {r['v']:.4f} launch {r['vf']:.4f} | clock {r['clock_range'][0]:+.2f}..{r['clock_range'][1]:+.2f} norm {r['norm_max']:.0e}")
    print(f"A_tiled {100*At:.1f}%  R = {Rr:.2f} (boot 68%: {np.percentile(bR,16):.2f}-{np.percentile(bR,84):.2f}; NaN-free {np.isfinite(bR).mean():.2f}) -> {verdict}")
    print(f"S: s5 mean / tiled mean = {Sratio:.3f} ({'PASS' if abs(Sratio-1) <= 0.10 else 'FAIL'})")
    dev = {net: np.array(row[net]['m']) / row[net]['mean'] - 1 for net in row}
    print('pattern correlations: s5-s3 %.2f, s5-s4 %.2f, s3-s4 %.2f' % (np.corrcoef(dev['glass5'], dev['glass'])[0, 1], np.corrcoef(dev['glass5'], dev['glass4'])[0, 1], np.corrcoef(dev['glass'], dev['glass4'])[0, 1]))
    res[point] = dict(R=float(Rr), R_boot68=[float(np.percentile(bR, 16)), float(np.percentile(bR, 84))], verdict=verdict, A_tiled=float(At), S_ratio=float(Sratio),
                      nets={k: {kk: vv for kk, vv in v.items() if kk != 'boot'} | dict(A_boot68=[float(np.percentile(v['boot'], 16)), float(np.percentile(v['boot'], 84))]) for k, v in row.items()})
json.dump(res, open('score_4k.json', 'w'), indent=1)
