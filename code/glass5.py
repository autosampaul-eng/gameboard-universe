"""Exp 4k: build glass s5 (13,824,000 nodes, ncell 120, seed 5) in stages, then measure c.
python glass5.py build|relax|finish --state DIR --deadline UNIXTIME"""
import numpy as np, os, sys, json, time, argparse
import crn, bigglass as B
ap = argparse.ArgumentParser(); ap.add_argument('mode'); ap.add_argument('--state', default='../state'); ap.add_argument('--deadline', type=float, default=1e30); a = ap.parse_args()
ncell, seed = 120, 5; aa = 4 * crn.D / np.sqrt(3); Lb = ncell * aa; n = 8 * ncell ** 3; tag = f'glass_n{n}_s{seed}'
NET = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'networks'); S = a.state; os.makedirs(S, exist_ok=True)
log = lambda s: print(f"[{time.strftime('%H:%M:%S')}] {s}", flush=True)
def save(pos, nb, v, st):
    np.save(f'{S}/pos.npy', pos); np.save(f'{S}/nb.npy', nb.astype(np.int32)); np.save(f'{S}/v.npy', v); np.save(f'{S}/st.npy', st)
def load():
    return np.load(f'{S}/pos.npy'), np.load(f'{S}/nb.npy').astype(np.int64), np.load(f'{S}/v.npy'), np.load(f'{S}/st.npy')
if a.mode == 'build':
    t0 = time.time(); pos, nb = B.random_network_lean(n, Lb, seed, passes=400); log(f'network built in {time.time()-t0:.0f}s')
    v = np.zeros_like(pos); st = B.new_state()
    B.relax(pos, nb, Lb, v, st, budget_s=a.deadline - time.time(), log=log); save(pos, nb, v, st); log(f'saved state, st={st.tolist()}')
elif a.mode == 'relax':
    pos, nb, v, st = load(); log(f'resuming at FIRE it {int(st[3])}, converged flag {st[6]}')
    if st[6] == 0: B.relax(pos, nb, Lb, v, st, budget_s=a.deadline - time.time(), log=log); save(pos, nb, v, st)
    log(f'st={st.tolist()}')
elif a.mode == 'finish':
    pos, nb, v, st = load(); assert st[6] != 0, 'relaxation not finished'
    from numba import njit
    @njit
    def forces(pos, nb, Lb):
        g = np.zeros_like(pos)
        for i in range(pos.shape[0]): crn.atom_eg(i, pos, nb, Lb, g)
        return np.abs(g).max(), np.sqrt((g ** 2).sum(1)).mean()
    fm, fmean = forces(pos, nb, Lb)
    q = crn.quality(pos, nb, Lb); q.update(Lb=Lb, fire_iterations=int(st[3]), converged=bool(st[6] == 1), fmax=float(fm), fmean=float(fmean))
    log(json.dumps(q)); json.dump(q, open(f'{NET}/{tag}_raw.json', 'w'))
    np.save(f'{NET}/{tag}_raw_pos.npy', pos); np.save(f'{NET}/{tag}_raw_nb.npy', nb.astype(np.int32)); del v
    from exp4 import load_net, adj, omega
    eta = 0.05; out = {}
    pos, nb, rev, d, LB = load_net('glass5'); A = adj(nb); ce = np.cos(eta)
    nC = 1500; j = np.arange(nC); xs = np.cos(np.pi * (j + 0.5) / nC); fv = np.sqrt(1 - ce ** 2 * xs ** 2)
    coef = np.array([2 / nC * (fv * np.cos(qq * np.pi * (j + 0.5) / nC)).sum() for qq in range(nC)]); coef[0] /= 2
    nuse = int(np.nonzero(np.abs(coef) > 1e-10)[0].max()) + 1
    g = 2 * np.pi / LB; c0 = np.array([LB / 2] * 3); r = pos - c0; r -= LB * np.round(r / LB)
    for dname, kv in [('100n1', g * np.array([1., 0, 0])), ('100n2', g * np.array([2., 0, 0])), ('100n3', g * np.array([3., 0, 0])), ('100n4', g * np.array([4., 0, 0])), ('123m2', 2 * g * np.array([1., 2, 3]))]:
        kk = np.linalg.norm(kv); ph0 = np.exp(-1j * (r @ kv))
        S0 = np.exp(-(r ** 2).sum(1) / (4 * 120.0 ** 2)) / ph0
        Bf = lambda v: (A @ v) / 4; T0, T1 = S0, Bf(S0); acc = coef[0] * T0 + coef[1] * T1
        for qq in range(2, nuse): T0, T1 = T1, 2 * Bf(T1) - T0; acc = acc + coef[qq] * T1
        del T0, T1; Sp, Sc = S0, ce * Bf(S0) + 1j * acc; del acc
        ph = [np.angle((Sp * ph0).sum()), np.angle((Sc * ph0).sum())]
        for t in range(2, 81):
            Sn = (ce / 2) * (A @ Sc) - Sp; Sp, Sc = Sc, Sn; ph.append(np.angle((Sc * ph0).sum()))
        ph = np.unwrap(ph); w = np.polyfit(np.arange(len(ph)), ph, 1)[0]
        c_true = np.sqrt(max(w ** 2 - eta ** 2, 0)) / kk; wmf = omega(d, kv, eta); c_mf = np.sqrt(wmf ** 2 - eta ** 2) / kk
        out[dname] = dict(k=float(kk), w=float(w), c_true=float(c_true), c_meanfield=float(c_mf)); log(f'{dname} k {kk:.5f} c_true {c_true:.5f} c_meanfield {c_mf:.5f}')
    out['c'] = out['100n2']['c_true']; out['direction_spread'] = abs(out['123m2']['c_true'] / out['100n2']['c_true'] - 1)
    json.dump(out, open(f'{NET}/{tag}_c.json', 'w'), indent=1); log(f"c = {out['c']:.5f} spread {out['direction_spread']:.4f}")
