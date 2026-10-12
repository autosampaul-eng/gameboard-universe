import json, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
S = json.load(open('score_4k.json')); DIRS = ['100', '010', '001', '123', '312', '231']; lab = ['[100]', '[010]', '[001]', '[1,2,3]', '[3,1,2]', '[2,3,1]']
nets = [('glass', 's3, tiled 2×2×2', '#20808D', 'o', -0.18), ('glass4', 's4, tiled 2×2×2', '#A84B2F', 's', 0.0), ('glass5', 's5, untiled (8× independent volume)', '#1B474D', 'D', 0.18)]
fig, axs = plt.subplots(1, 2, figsize=(13.5, 5.0), sharey=True)
for ax, pt, title in [(axs[0], '5', 'Long wavelength (35.6 bonds): spread shrinks in the untiled glass'), (axs[1], '1', 'Short wavelength (17.8 bonds): all spreads small')]:
    for net, name, col, m, off in nets:
        r = S[pt]['nets'][net]; dev = 100 * (np.array(r['m']) / r['mean'] - 1)
        raw = np.std(dev, ddof=1)
        ax.plot(np.arange(6) + off, dev, m, color=col, ms=8, label=f"{name}: spread {raw:.1f}%")
    ax.axhline(0, color='#7A7974', lw=0.8); ax.set_xticks(range(6)); ax.set_xticklabels(lab)
    ax.set_title(title, fontsize=10.5); ax.spines[['top', 'right']].set_visible(False); ax.legend(frameon=False, fontsize=9, loc='upper right')
    ax.set_xlabel('direction of travel')
axs[0].set_ylabel('scattering rate vs the glass\'s own average (%)'); axs[0].set_ylim(-13, 16)
plt.tight_layout(); plt.savefig('exp4k_results.png', dpi=150)
