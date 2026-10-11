"""Build the replication glass s4 with exactly the recipe of make_raw_big.py (ncell 60, seed 4), then measure its speed limit c
by the Exp 4c method (diag_c3): node-total wave equation, eta = 0.05, longest box-fitting waves."""
import numpy as np, time, json, os, sys
import crn, build
NET = os.environ.get('GBU_NETWORKS', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'networks'))
ncell, seed = 60, 4
a = 4 * crn.D / np.sqrt(3); Lb = ncell * a; n = 8 * ncell ** 3; tag = f'glass_n{n}_s{seed}'
T0 = time.time(); pos, nb = build.random_network(n, Lb, seed); print('built', time.time() - T0, flush=True)
crn.global_relax(pos, nb, Lb, 1e-3, 5000)
q = crn.quality(pos, nb, Lb); q['secs'] = time.time() - T0; q['Lb'] = Lb; print(json.dumps(q), flush=True)
np.save(f'{NET}/{tag}_raw_pos.npy', pos); np.save(f'{NET}/{tag}_raw_nb.npy', nb.astype(np.int32)); json.dump(q, open(f'{NET}/{tag}_raw.json', 'w'))

from exp4 import load_net, adj, mu_of_k, omega
eta = 0.05; out = {}
pos, nb, rev, d, LB = load_net('glass4'); A = adj(nb); ce = np.cos(eta)
nC = 1500; j = np.arange(nC); xs = np.cos(np.pi * (j + 0.5) / nC); fv = np.sqrt(1 - ce ** 2 * xs ** 2)
coef = np.array([2 / nC * (fv * np.cos(qq * np.pi * (j + 0.5) / nC)).sum() for qq in range(nC)]); coef[0] /= 2
nuse = int(np.nonzero(np.abs(coef) > 1e-10)[0].max()) + 1
g = 2 * np.pi / LB
for dname, kv in [('100n1', g * np.array([1., 0, 0])), ('100n2', g * np.array([2., 0, 0])), ('100n3', g * np.array([3., 0, 0])), ('123m1', g * np.array([1., 2, 3]))]:
    kk = np.linalg.norm(kv); c0 = np.array([LB / 2] * 3); r = pos - c0; r -= LB * np.round(r / LB)
    S0 = np.exp(-(r ** 2).sum(1) / (4 * 60.0 ** 2)) * np.exp(1j * (r @ kv))
    B = lambda v: (A @ v) / 4; T0_, T1 = S0, B(S0); acc = coef[0] * T0_ + coef[1] * T1
    for qq in range(2, nuse): T0_, T1 = T1, 2 * B(T1) - T0_; acc = acc + coef[qq] * T1
    Sp, Sc = S0, ce * B(S0) + 1j * acc
    ph = [np.angle((Sp * np.exp(-1j * (r @ kv))).sum()), np.angle((Sc * np.exp(-1j * (r @ kv))).sum())]
    for t in range(2, 81):
        Sn = (ce / 2) * (A @ Sc) - Sp; Sp, Sc = Sc, Sn; ph.append(np.angle((Sc * np.exp(-1j * (r @ kv))).sum()))
    ph = np.unwrap(ph); w = np.polyfit(np.arange(len(ph)), ph, 1)[0]
    c_true = np.sqrt(max(w ** 2 - eta ** 2, 0)) / kk; wmf = omega(d, kv, eta); c_mf = np.sqrt(wmf ** 2 - eta ** 2) / kk
    out[dname] = dict(k=float(kk), w=float(w), c_true=float(c_true), c_meanfield=float(c_mf))
    print(dname, 'c_true %.5f c_meanfield %.5f ratio %.4f' % (c_true, c_mf, c_true / c_mf), flush=True)
out['c'] = out['100n1']['c_true']; out['direction_spread'] = abs(out['123m1']['c_true'] / out['100n1']['c_true'] - 1)
json.dump(out, open(f'{NET}/{tag}_c.json', 'w'), indent=1); print('c =', out['c'], 'spread', out['direction_spread'])
