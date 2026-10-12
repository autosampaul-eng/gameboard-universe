"""Exp 4j scoring, as registered in DESIGN_exp4j.md. Usage: python score_4j.py RESULTS_DIR RES_4I_DIR RES_4I_REP_DIR"""
import json, glob, sys, numpy as np
from scipy.optimize import curve_fit
D, D4I, D4R = sys.argv[1:4]
L = lambda d: [json.load(open(f)) for f in sorted(glob.glob(f'{d}/job_*.json'))]
J = L(D); I3 = L(D4I); I4 = L(D4R)
BIGTEST_V = 0.5494640716399202; K_REF = 0.07675          # bigtest diamond x3 point 5; middle wavelength (j index 3)
out = {}
def ms(x): x = np.asarray(x, float); return float(x.mean()), float(x.std(ddof=1) / np.sqrt(len(x)))
# --- 0a regression
reg = [r for r in J if r['set'] == 'REG']; d24 = [r for r in J if r['set'] == 'C' and r['i'] == 2 and r['j'] == 4]
e1 = abs(reg[0]['coh']['v'] - BIGTEST_V) if reg else np.inf; e2 = abs(d24[0]['coh']['v'] - BIGTEST_V) if d24 else np.inf
out['0a'] = dict(reg_diff=e1, new_path_diff=e2, pass_=bool(e1 < 1e-10 and e2 < 1e-10))
print(f"0a regression: default path |dv| {e1:.1e}, new path (eta=0.15 passed) |dv| {e2:.1e} -> {'PASS' if out['0a']['pass_'] else 'FAIL'}")
# --- 0b diamond controls
bad = []
for r in [r for r in J if r['set'] == 'C']:
    ok = r['ds_max'] < 0.5 and abs(r['inv_tau']) < 2e-5 and abs(r['coh']['v'] / r['whole']['v'] - 1) < 1e-3
    if not ok: bad.append((r['i'], r['j']))
out['0b'] = dict(failed=bad, not_working=len(bad) >= 4)
print(f"0b diamond controls: {21 - len(bad)}/21 pass; failed cells {bad}" + ("  -> EXPERIMENT NOT WORKING" if len(bad) >= 4 else ''))
# --- cells
G = [r for r in J if r['set'] == 'G']; cells = {}
for i in range(3):
    for j in range(7):
        g = [r for r in G if r['i'] == i and r['j'] == j]
        if len(g) == 0: continue
        cells[(i, j)] = dict(n=len(g), k0=g[0]['k0'], gb=g[0]['gb'], vf=g[0]['vf'], delta=ms([r['delta'] for r in g]), inv_ell=ms([r['inv_ell'] for r in g]),
                             raw=ms([r['shortfall_raw'] for r in g]), delta_spec=ms([r['delta_spec'] for r in g]), v=ms([r['coh']['v'] for r in g]),
                             clock=[min(r['clock_err'] for r in g), max(r['clock_err'] for r in g)], norm=max(r['norm'] for r in g))
# --- 0c box check vs 4i / 4i-R at the eta = 0.15 cells
box_ok = True; out['0c'] = []
for i in range(3):
    j = 2 * i; pt = 2 * i + 1
    for net, old in (('glass', I3), ('glass4', I4)):
        new = [r for r in G if r['i'] == i and r['j'] == j and r['net'] == net]; o = [r for r in old if r.get('set') == 'G' and r['point'] == pt]
        if not new or not o: continue
        for name, f in (('1/ell', lambda r: r['inv_ell']), ('v_p', lambda r: r['coh']['v'])):
            a, sa = ms([f(r) for r in new]); b, sb = ms([f(r) for r in o]); z = abs(a - b) / np.hypot(sa, sb)
            box_ok &= bool(z < 2); out['0c'].append(dict(net=net, i=i, j=j, q=name, new=a, old=b, z=float(z)))
            print(f"0c box check {net} speed {i} wavelength {j+1} {name}: 3x3x3 {a:.5g} vs 2x2x2 {b:.5g} ({z:.1f} SE)")
print('0c:', 'PASS' if box_ok else 'FAIL (box size matters; Q1-Q3 labelled unreliable)')
use = [(i, j) for (i, j) in cells if (i, j) not in bad]
k = np.array([cells[c]['k0'] for c in use]); gb = np.array([cells[c]['gb'] for c in use])
# --- Q1 / Q3 (delta): delta = A k^q gb^s, weighted least squares on plain values
dl = np.array([cells[c]['delta'][0] for c in use]); de = np.array([cells[c]['delta'][1] for c in use])
f1 = lambda X, A, q, s: A * (X[0] / K_REF) ** q * X[1] ** s
p0 = [max(dl.mean(), 1e-3), 1.5, 0.0]
P, Cv = curve_fit(f1, np.vstack([k, gb]), dl, p0=p0, sigma=de, absolute_sigma=True, maxfev=20000); E = np.sqrt(np.diag(Cv))
chi1 = float((((dl - f1(np.vstack([k, gb]), *P)) / de) ** 2).sum() / (len(dl) - 3))
q, qe, s, se_ = P[1], E[1], P[2], E[2]
if chi1 >= 3: v1 = 'not a single power law'
elif abs(q - 2) <= 2 * qe: v1 = 'ordinary lattice dispersion'
elif q < 2: v1 = 'fades more slowly than an ordinary lattice'
else: v1 = 'fades faster than an ordinary lattice'
v3d = 'depends on speed' if abs(s) > 2 * se_ else 'no speed dependence detected'
print(f"\nQ1 excess shortfall: q = {q:.3f} ± {qe:.3f}, chi2/dof {chi1:.2f} -> {v1}  (A at k_ref, gb=1: {100*P[0]:.3f}%)")
print(f"Q3 (shortfall): s = {s:.3f} ± {se_:.3f} -> {v3d}")
# --- Q2 / Q3 (scattering): ln(1/ell) = a + p x + b x^2 + s' ln gb, x = ln(k/K_REF)
y = np.array([cells[c]['inv_ell'][0] for c in use]); ye = np.array([cells[c]['inv_ell'][1] for c in use])
x = np.log(k / K_REF); Y = np.log(y); W = (y / ye) ** 2
A_ = np.vstack([np.ones_like(x), x, x ** 2, np.log(gb)]).T; Cs = np.linalg.inv(A_.T @ (A_ * W[:, None])); B = Cs @ (A_.T @ (W * Y))
chi2 = float((W * (Y - A_ @ B) ** 2).sum() / (len(Y) - 4)); Bse = np.sqrt(np.diag(Cs))
xl = np.log(K_REF * 0 + cells[(0, 5)]['k0'] / K_REF) if (0, 5) in cells else np.log(0.05427 / K_REF)
slope_long = B[1] + 2 * B[2] * xl; slope_long_se = float(np.sqrt(np.array([0, 1, 2 * xl, 0]) @ Cs @ np.array([0, 1, 2 * xl, 0])))
v2 = 'steepens toward long wavelengths' if B[2] < -2 * Bse[2] else 'no bend detected'
v3s = 'depends on speed' if abs(B[3]) > 2 * Bse[3] else 'no speed dependence detected'
print(f"Q2 scattering: p (at 29.9 bonds) = {B[1]:.3f} ± {Bse[1]:.3f}, b = {B[2]:.3f} ± {Bse[2]:.3f}, chi2/dof {chi2:.2f} -> {v2}; "
      f"local slope at 42.3 bonds {slope_long:.2f} ± {slope_long_se:.2f} (4 within 2 SE: {abs(slope_long-4) <= 2*slope_long_se})")
print(f"Q3 (scattering): s' = {B[3]:.3f} ± {Bse[3]:.3f} -> {v3s}")
# --- clocks
cl = [r['clock_err'] for r in G]; nm = [r['norm'] for r in G]
print(f"Clocks {min(cl):+.2f}..{max(cl):+.2f}% ({'PASS' if max(abs(c) for c in cl) < 3 else 'FAIL'}); total {max(nm):.1e} ({'PASS' if max(nm) < 1e-9 else 'FAIL'})")
# --- reported, not scored
print('\nPer cell (glass mean over 6 runs):')
for (i, j), c in sorted(cells.items()):
    print(f"  {c['vf']:.3f}c j{j+1} k0 {c['k0']:.5f}: delta {100*c['delta'][0]:+.3f} ± {100*c['delta'][1]:.3f}% (raw {100*c['raw'][0]:+.3f}%, spectral {100*c['delta_spec'][0]:+.3f}%) "
          f"1/ell {c['inv_ell'][0]:.4e} ± {c['inv_ell'][1]:.1e} | clock {c['clock'][0]:+.2f}..{c['clock'][1]:+.2f} norm {c['norm']:.0e}")
for i in range(3):
    m = [n for n, c in enumerate(use) if c[0] == i]
    if len(m) >= 3:
        try:
            Pi, Ci = curve_fit(lambda kk, A, qq: A * (kk / K_REF) ** qq, k[m], dl[m], p0=[max(dl[m].mean(), 1e-3), 1.5], sigma=de[m], absolute_sigma=True, maxfev=20000)
            pi = np.polyfit(np.log(k[m]), np.log(y[m]), 1, w=y[m] / ye[m]); print(f"  speed {i}: q = {Pi[1]:.2f} ± {np.sqrt(Ci[1,1]):.2f}; scattering slope p = {pi[0]:.2f}")
        except Exception as ex: print('  speed', i, 'fit failed', ex)
try:
    Pf, Cf = curve_fit(lambda X, A, qq, ss, Bf: A * (X[0] / K_REF) ** qq * X[1] ** ss + Bf, np.vstack([k, gb]), dl, p0=[*P, 0.0], sigma=de, absolute_sigma=True, maxfev=40000)
    print(f"  floor fit: q {Pf[1]:.2f} ± {np.sqrt(Cf[1,1]):.2f}, floor B = {100*Pf[3]:+.3f} ± {100*np.sqrt(Cf[3,3]):.3f}%")
except Exception as ex: print('  floor fit failed', ex)
dm = [r for r in J if r['set'] == 'C' and (r['i'], r['j']) not in bad]
try:
    kd = np.array([r['k0'] for r in dm]); gd = np.array([r['gb'] for r in dm]); dd = np.array([r['delta'] for r in dm])
    Pd, Cd = curve_fit(f1, np.vstack([kd, gd]), dd, p0=[max(dd.mean(), 1e-3), 2, 0], maxfev=20000)
    print(f"  diamond (reference, unweighted single runs): q {Pd[1]:.2f}, s {Pd[2]:.2f}, delta range {100*dd.min():+.3f}..{100*dd.max():+.3f}%")
except Exception as ex: print('  diamond fit failed', ex)
out.update(Q1=dict(q=float(q), q_se=float(qe), s=float(s), s_se=float(se_), A=float(P[0]), chi2dof=chi1, verdict=v1, verdict_speed=v3d),
           Q2=dict(p=float(B[1]), p_se=float(Bse[1]), b=float(B[2]), b_se=float(Bse[2]), s=float(B[3]), s_se=float(Bse[3]), chi2dof=chi2, slope_long=float(slope_long),
                   slope_long_se=slope_long_se, verdict=v2, verdict_speed=v3s, coef=B.tolist()),
           cells={f'{i},{j}': c for (i, j), c in cells.items()})
json.dump(out, open('score_4j.json', 'w'), indent=1, default=lambda o: o if not isinstance(o, np.generic) else o.item())
