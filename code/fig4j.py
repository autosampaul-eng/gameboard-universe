import json, glob, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
S = json.load(open('score_4j.json')); J = [json.load(open(f)) for f in sorted(glob.glob('res_4j/job_*.json'))]
K = 0.07675; cols = ['#20808D', '#A84B2F', '#1B474D']; mk = ['o', 's', 'D']; spd = ['0.800c', '0.686c', '0.555c']
fig, (a1, a2) = plt.subplots(1, 2, figsize=(13.5, 5.2))
Q1 = S['Q1']; Q2 = S['Q2']
for i in range(3):
    c = [S['cells'][f'{i},{j}'] for j in range(7)]; k = np.array([x['k0'] for x in c]); lam = 2 * np.pi / k / 2.735; gb = c[0]['gb']
    d = 100 * np.array([x['delta'][0] for x in c]); de = 100 * np.array([x['delta'][1] for x in c])
    a1.errorbar(lam, d, yerr=de, fmt=mk[i], color=cols[i], ms=7, capsize=3, label=f'glass, {spd[i]}')
    kk = np.linspace(k.min() * 0.95, k.max() * 1.05, 50); a1.plot(2 * np.pi / kk / 2.735, 100 * Q1['A'] * (kk / K) ** Q1['q'] * gb ** Q1['s'], '-', color=cols[i], lw=1.1, alpha=0.8)
    y = np.array([x['inv_ell'][0] for x in c]); ye = np.array([x['inv_ell'][1] for x in c])
    a2.errorbar(lam * (1 + 0.012 * (i - 1)), 1 / y / 2.735, yerr=(ye / y) / y / 2.735, fmt=mk[i], color=cols[i], ms=7, capsize=3, label=f'glass, {spd[i]}')
dm = sorted([r for r in J if r['set'] == 'C' and r['i'] == 2], key=lambda r: r['j'])
a1.plot([2 * np.pi / r['k0'] / 2.735 for r in dm], [100 * r['delta'] for r in dm], 'D', mfc='none', color='#7A7974', ms=7, label='diamond crystal, 0.555c')
lam_ref = 2 * np.pi / np.array([0.12907, 0.04563]) / 2.735; ref = 100 * Q1['A'] * (np.array([0.12907, 0.04563]) / K) ** 2 * (4 / 3) ** Q1['s']
a1.plot(lam_ref, ref, ':', color='#7A7974', lw=1.2, label='slope of an ordinary lattice (1/λ²)')
x = np.linspace(np.log(0.04563 / K) - 0.05, np.log(0.12907 / K) + 0.05, 50); B = Q2['coef']
a2.plot(2 * np.pi / (K * np.exp(x)) / 2.735, 1 / np.exp(B[0] + B[1] * x + B[2] * x ** 2 + B[3] * np.log(0.943)) / 2.735, '-', color='#7A7974', lw=1.1, label=f"fit: 1/ℓ ∝ k^{Q2['p']:.2f}, no speed term needed")
for a in (a1, a2):
    a.set_xscale('log'); a.set_yscale('log'); a.set_xticks([18, 21, 25, 30, 36, 42, 50]); a.set_xticklabels(['18', '21', '25', '30', '36', '42', '50']); a.minorticks_off()
    a.set_xlabel('wavelength (bonds)'); a.spines[['top', 'right']].set_visible(False)
a1.set_yticks([0.1, 0.2, 0.5, 1, 2, 5, 10]); a1.set_yticklabels(['0.1', '0.2', '0.5', '1', '2', '5', '10']); a1.set_ylim(0.08, 14)
a1.set_ylabel('excess speed shortfall δ (%)'); a1.set_title(f"Extra slowdown fades as 1/λ^{Q1['q']:.1f}; slower particles are slowed more", fontsize=10.5)
a1.legend(frameon=False, fontsize=8.5, loc='lower left')
a2.set_yticks([500, 1000, 2000, 5000, 10000, 20000]); a2.set_yticklabels(['500', '1,000', '2,000', '5,000', '10,000', '20,000'])
a2.set_ylabel('mean free path (bonds)'); a2.set_title('Scattering depends on wavelength only, not speed', fontsize=10.5); a2.legend(frameon=False, fontsize=8.5, loc='upper left')
plt.tight_layout(); plt.savefig('exp4j_results.png', dpi=150)
